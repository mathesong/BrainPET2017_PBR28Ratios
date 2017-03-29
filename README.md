``` r
library(tidyverse)
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(stringr)
library(kinfitr)
```

Read in the demographic data
----------------------------

``` r
demog <- readxl::read_excel('../RawData/TrT_chemistry_demograph.xlsx') %>%
  select(Subjname=Akronym, Gender=Sex, Age, Genotype, 
         PET_same_day, `MBq PET1`, `MBq PET2`, bodyMass=Weight_kg) %>%
  tidyr::gather(PETNo, injRad, dplyr::contains('MBq')) %>%
  mutate(PETNo = as.numeric(stringr::str_match(PETNo, '\\d'))) %>%
  mutate(PET = paste(Subjname, PETNo, sep='_'))
```

``` r
demog %>%
  select(Age, InjectedRadioactivity = injRad) %>%
  psych::describe() %>%
  pander::pandoc.table(digits=3, caption = "Summary Statistics", split.tables=Inf)
```

    ## 
    ## -------------------------------------------------------------------------------------------------------------------------------
    ##           &nbsp;             vars   n    mean    sd    median   trimmed   mad    min   max   range    skew    kurtosis    se   
    ## --------------------------- ------ ---- ------ ------ -------- --------- ------ ----- ----- ------- -------- ---------- -------
    ##           **Age**             1     24   23.9   2.99     24      23.8     3.71   20    29      9     0.226     -1.32     0.611 
    ## 
    ##  **InjectedRadioactivity**    2     24   395    51.7    403       401     51.9   248   462    214    -0.997    0.722     10.6  
    ## -------------------------------------------------------------------------------------------------------------------------------
    ## 
    ## Table: Summary Statistics

``` r
counts <- demog %>%
  filter(PETNo==1)

pander::pandoc.table(table(counts$Gender, counts$Genotype), caption = "Gender and Genotype")
```

    ## 
    ## --------------------
    ##  &nbsp;   HAB   MAB 
    ## -------- ----- -----
    ##  **F**     2     4  
    ## 
    ##  **M**     4     2  
    ## --------------------
    ## 
    ## Table: Gender and Genotype

Read in the TACs and the blood data
-----------------------------------

``` r
dat <- readRDS('../RawData/pbrdat.rds')

datdf <- map(dat, 'tacdf') %>%
  bind_rows(.id = "id") %>%
  select(PET = id, Times = Times, Weights=weights, FC = FSLSFC, TC = FSLSTC, 
         STR = FSLSSTR, THA = FSLSTHA, WB=WB, CBL=FSLSCER) %>%
  group_by(PET) %>%
  nest() %>%
  rename(tacs=data) %>%
  mutate(Subjname = stringr::str_extract(PET, "(^[a-z]*)"), 
         PETNo = as.numeric(stringr::str_extract(PET, "\\d$")),
         input=map(dat, 'input')) %>%
  inner_join(demog)
```

    ## Joining, by = c("PET", "Subjname", "PETNo")

Fit the delay
-------------

``` r
datdf <- datdf %>%
  group_by(PET) %>%
  mutate(WB_delay = map2(tacs, input, ~twotcm(t_tac = .x$Times/60, tac = .x$WB, 
                                              input = .y, frameStartEnd = c(1,33), inpshift.upper = 1)))

saveRDS(datdf, file='datdf.rds')
```

Convert the data to longer format
---------------------------------

``` r
tacs <- datdf %>%
  select(PET, tacs, Subjname, PETNo) %>%
  unnest() %>%
  gather(Region, TAC, -PET, -Subjname, -PETNo, -Weights, -Times) %>%
  group_by(PET, Subjname, PETNo, Region) %>%
  nest() %>%
  rename(tacs=data)
  

longdat <- datdf %>%
  select(PET, Subjname, PETNo, input, WB_delay, bodyMass, injRad) %>%
  inner_join(tacs)
```

    ## Joining, by = c("PET", "Subjname", "PETNo")

Define functions for fitting the models
---------------------------------------

``` r
# Total SUV
calcSUV <- function(tacs, bodymass, injRad, frameStartEnd) {
  SUV(t_tac = tacs$Times/60, 
         tac = tacs$TAC*0.037,    # To kBq - because of kg
         injRad = injRad,
         bodymass = bodymass,
         frameStartEnd=c( frameStartEnd[1],frameStartEnd[2]))
}

# 40-60 Minute SUV
calcSUV_4060 <- function(SUVout) {
  
    interptime = seq(SUVout$tacs$Time[1], rev(SUVout$tacs$Time)[1], by=1/60)
    interptac  = pracma::interp1(x = SUVout$tacs$Time, SUVout$tacs$SUV, xi = interptime, method = 'linear')
    step=interptime[2] - interptime[1]
    
    SUVdf <- data.frame(interptime, interptac) %>%
      filter(interptime > 40 & interptime <= 60)
    
    out <- mean(SUVdf$interptac)
    return(out)
}

# MA1 using the fitted delay and vB from delayFit
fitma1 <- function(tacs, input, delayFit) {
  ma1(t_tac = tacs$Times/60, tac = tacs$TAC, input = input, tstarIncludedFrames = 6,
      inpshift = delayFit$par$inpshift, frameStartEnd=c(1,33), weights=tacs$Weights,
      vB=delayFit$par$vB)
}

# 2TCM using the fitted delay and vB from delayFit
fit2tcm <- function(tacs, input, delayFit) {
  twotcm(t_tac = tacs$Times/60, tac = tacs$TAC, input = input, inpshift = delayFit$par$inpshift,
         vB=delayFit$par$vB, frameStartEnd=c(1,33), weights=tacs$Weights)
}
```

Fit all the models
------------------

``` r
longdat <- longdat %>%
  # SUV Total
  mutate(suvout_tot = pmap(list(tacs, bodyMass, injRad), calcSUV, frameStartEnd=c(1,33))) %>%
  mutate(SUV_tot = purrr::map(suvout_tot, c('par', 'intSUV')) %>% do.call(rbind, .) %>% as.numeric()) %>%
  
  # SUV 40-60 Minutes
  mutate(SUV4060 = purrr::map(suvout_tot, calcSUV_4060) %>% do.call(rbind, .) %>% as.numeric()) %>%
  
  # MA1
  mutate(fit_ma1 = pmap(list(tacs, input, WB_delay), fitma1)) %>%
  mutate(Vt_ma1 = purrr::map(fit_ma1, c('par', 'Vt')) %>% do.call(rbind, .) %>% as.numeric()) %>%
    
  # 2TCM using fitted vB and delay
  mutate(fit_2tcm= pmap(list(tacs, input, WB_delay), fit2tcm)) %>%
  mutate(Vt_2tcm = purrr::map(fit_2tcm, c('par', 'Vt')) %>% do.call(rbind, .) %>% as.numeric())
  

saveRDS(longdat, file='longdat.rds')
```

Calculate the ratio metrics for all measures
--------------------------------------------

``` r
SUV_tot <- longdat %>%
  ungroup() %>%
  select(Subjname, PETNo, Region, SUV_tot) %>%
  spread(Region, SUV_tot) %>%
  mutate(FC_CBL = FC/CBL, FC_WB=FC/WB) %>%
  mutate(Measure='SUV_Tot')

SUV_4060 <- longdat %>%
  ungroup() %>%
  select(Subjname, PETNo, Region, SUV4060) %>%
  spread(Region, SUV4060) %>%
  mutate(FC_CBL = FC/CBL, FC_WB=FC/WB) %>%
  mutate(Measure='SUV_4060')

VT_2tcm <- longdat %>%
  ungroup() %>%
  select(Subjname, PETNo, Region, Vt_2tcm) %>%
  spread(Region, Vt_2tcm) %>%
  mutate(FC_CBL = FC/CBL, FC_WB=FC/WB) %>%
  mutate(Measure='Vt_2tcm')

VT_MA1 <- longdat %>%
  ungroup() %>%
  select(Subjname, PETNo, Region, Vt_ma1) %>%
  spread(Region, Vt_ma1) %>%
  mutate(FC_CBL = FC/CBL, FC_WB=FC/WB) %>%
  mutate(Measure='Vt_ma1')

trtdata <- bind_rows(list(SUV_tot, SUV_4060, VT_2tcm, VT_MA1))
```

Calculate the Test-Retest Metrics
---------------------------------

### HABs and MABs together

``` r
trtout <- trtdata %>%
  gather(Region, Binding, -Subjname, -PETNo, -Measure) %>%
  filter(Subjname != 'mahi') %>%
  spread(PETNo, Binding) %>%
  group_by(Region, Measure) %>%
  do(trt = granviller::trt(.$`1`, .$`2`)$tidy) %>%
  unnest() %>%
  arrange(Measure, Region)

pander::pandoc.table(trtout, digits=2, split.tables=Inf, caption='Test-Retest Analysis for HABs and MABs together')
```

<table style="width:100%;">
<caption>Test-Retest Analysis for HABs and MABs together</caption>
<colgroup>
<col width="7%" />
<col width="9%" />
<col width="6%" />
<col width="7%" />
<col width="7%" />
<col width="7%" />
<col width="7%" />
<col width="9%" />
<col width="6%" />
<col width="6%" />
<col width="7%" />
<col width="5%" />
<col width="13%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Region</th>
<th align="center">Measure</th>
<th align="center">mean</th>
<th align="center">sd</th>
<th align="center">cov</th>
<th align="center">se</th>
<th align="center">skew</th>
<th align="center">kurtosis</th>
<th align="center">icc</th>
<th align="center">aapd</th>
<th align="center">sem</th>
<th align="center">md</th>
<th align="center">avgpercchange</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">SUV_4060</td>
<td align="center">0.98</td>
<td align="center">0.3</td>
<td align="center">0.3</td>
<td align="center">0.063</td>
<td align="center">0.5</td>
<td align="center">-1</td>
<td align="center">0.9</td>
<td align="center">10</td>
<td align="center">0.092</td>
<td align="center">26</td>
<td align="center">-3.3</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">SUV_4060</td>
<td align="center">0.92</td>
<td align="center">0.27</td>
<td align="center">0.3</td>
<td align="center">0.058</td>
<td align="center">0.57</td>
<td align="center">-0.92</td>
<td align="center">0.88</td>
<td align="center">13</td>
<td align="center">0.097</td>
<td align="center">29</td>
<td align="center">-4.3</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">SUV_4060</td>
<td align="center">0.95</td>
<td align="center">0.059</td>
<td align="center">0.062</td>
<td align="center">0.013</td>
<td align="center">-0.72</td>
<td align="center">0.24</td>
<td align="center">0.63</td>
<td align="center">4.4</td>
<td align="center">0.036</td>
<td align="center">10</td>
<td align="center">-1</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">SUV_4060</td>
<td align="center">1</td>
<td align="center">0.04</td>
<td align="center">0.039</td>
<td align="center">0.0085</td>
<td align="center">-0.24</td>
<td align="center">-1.2</td>
<td align="center">0.89</td>
<td align="center">1.5</td>
<td align="center">0.013</td>
<td align="center">3.5</td>
<td align="center">0.25</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">SUV_4060</td>
<td align="center">0.89</td>
<td align="center">0.27</td>
<td align="center">0.31</td>
<td align="center">0.058</td>
<td align="center">0.43</td>
<td align="center">-0.92</td>
<td align="center">0.84</td>
<td align="center">16</td>
<td align="center">0.11</td>
<td align="center">33</td>
<td align="center">-5.7</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">SUV_4060</td>
<td align="center">0.93</td>
<td align="center">0.27</td>
<td align="center">0.29</td>
<td align="center">0.057</td>
<td align="center">0.48</td>
<td align="center">-0.85</td>
<td align="center">0.89</td>
<td align="center">11</td>
<td align="center">0.089</td>
<td align="center">27</td>
<td align="center">-4.4</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">SUV_4060</td>
<td align="center">1.2</td>
<td align="center">0.37</td>
<td align="center">0.32</td>
<td align="center">0.08</td>
<td align="center">0.45</td>
<td align="center">-0.84</td>
<td align="center">0.87</td>
<td align="center">15</td>
<td align="center">0.13</td>
<td align="center">31</td>
<td align="center">-4.8</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">SUV_4060</td>
<td align="center">0.9</td>
<td align="center">0.24</td>
<td align="center">0.27</td>
<td align="center">0.052</td>
<td align="center">0.55</td>
<td align="center">-0.74</td>
<td align="center">0.86</td>
<td align="center">13</td>
<td align="center">0.092</td>
<td align="center">28</td>
<td align="center">-4.6</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">SUV_Tot</td>
<td align="center">80</td>
<td align="center">17</td>
<td align="center">0.21</td>
<td align="center">3.6</td>
<td align="center">0.53</td>
<td align="center">-0.86</td>
<td align="center">0.83</td>
<td align="center">11</td>
<td align="center">7.1</td>
<td align="center">24</td>
<td align="center">-5.6</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">SUV_Tot</td>
<td align="center">78</td>
<td align="center">17</td>
<td align="center">0.21</td>
<td align="center">3.5</td>
<td align="center">0.44</td>
<td align="center">-0.69</td>
<td align="center">0.8</td>
<td align="center">12</td>
<td align="center">7.5</td>
<td align="center">26</td>
<td align="center">-6.2</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">SUV_Tot</td>
<td align="center">0.98</td>
<td align="center">0.045</td>
<td align="center">0.046</td>
<td align="center">0.0097</td>
<td align="center">-0.29</td>
<td align="center">-1</td>
<td align="center">0.77</td>
<td align="center">2.5</td>
<td align="center">0.022</td>
<td align="center">6.1</td>
<td align="center">-0.58</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">SUV_Tot</td>
<td align="center">1.1</td>
<td align="center">0.039</td>
<td align="center">0.036</td>
<td align="center">0.0082</td>
<td align="center">-0.3</td>
<td align="center">-1.2</td>
<td align="center">0.9</td>
<td align="center">1.3</td>
<td align="center">0.012</td>
<td align="center">3.1</td>
<td align="center">0.079</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">SUV_Tot</td>
<td align="center">76</td>
<td align="center">17</td>
<td align="center">0.22</td>
<td align="center">3.6</td>
<td align="center">0.23</td>
<td align="center">-0.79</td>
<td align="center">0.77</td>
<td align="center">14</td>
<td align="center">8.1</td>
<td align="center">29</td>
<td align="center">-6.3</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">SUV_Tot</td>
<td align="center">75</td>
<td align="center">15</td>
<td align="center">0.2</td>
<td align="center">3.3</td>
<td align="center">0.5</td>
<td align="center">-0.53</td>
<td align="center">0.79</td>
<td align="center">12</td>
<td align="center">7</td>
<td align="center">26</td>
<td align="center">-6.3</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">SUV_Tot</td>
<td align="center">92</td>
<td align="center">20</td>
<td align="center">0.21</td>
<td align="center">4.2</td>
<td align="center">0.43</td>
<td align="center">-0.36</td>
<td align="center">0.74</td>
<td align="center">14</td>
<td align="center">10</td>
<td align="center">30</td>
<td align="center">-6.5</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">SUV_Tot</td>
<td align="center">72</td>
<td align="center">14</td>
<td align="center">0.19</td>
<td align="center">3</td>
<td align="center">0.45</td>
<td align="center">-0.44</td>
<td align="center">0.75</td>
<td align="center">13</td>
<td align="center">7</td>
<td align="center">27</td>
<td align="center">-6.2</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">Vt_2tcm</td>
<td align="center">3</td>
<td align="center">1.7</td>
<td align="center">0.55</td>
<td align="center">0.35</td>
<td align="center">1.2</td>
<td align="center">0.61</td>
<td align="center">0.93</td>
<td align="center">21</td>
<td align="center">0.45</td>
<td align="center">42</td>
<td align="center">16</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.9</td>
<td align="center">1.5</td>
<td align="center">0.53</td>
<td align="center">0.33</td>
<td align="center">1</td>
<td align="center">0.26</td>
<td align="center">0.93</td>
<td align="center">19</td>
<td align="center">0.42</td>
<td align="center">40</td>
<td align="center">13</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">Vt_2tcm</td>
<td align="center">0.98</td>
<td align="center">0.07</td>
<td align="center">0.071</td>
<td align="center">0.015</td>
<td align="center">-0.15</td>
<td align="center">-0.57</td>
<td align="center">0.54</td>
<td align="center">4.7</td>
<td align="center">0.048</td>
<td align="center">13</td>
<td align="center">-2.2</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">Vt_2tcm</td>
<td align="center">1</td>
<td align="center">0.037</td>
<td align="center">0.036</td>
<td align="center">0.0079</td>
<td align="center">-0.13</td>
<td align="center">-1.3</td>
<td align="center">0.52</td>
<td align="center">3</td>
<td align="center">0.026</td>
<td align="center">6.9</td>
<td align="center">-0.98</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.8</td>
<td align="center">1.5</td>
<td align="center">0.53</td>
<td align="center">0.31</td>
<td align="center">1.3</td>
<td align="center">1.1</td>
<td align="center">0.93</td>
<td align="center">19</td>
<td align="center">0.39</td>
<td align="center">40</td>
<td align="center">13</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.9</td>
<td align="center">1.6</td>
<td align="center">0.54</td>
<td align="center">0.34</td>
<td align="center">1.2</td>
<td align="center">0.69</td>
<td align="center">0.94</td>
<td align="center">19</td>
<td align="center">0.39</td>
<td align="center">37</td>
<td align="center">15</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">Vt_2tcm</td>
<td align="center">3.8</td>
<td align="center">2.3</td>
<td align="center">0.6</td>
<td align="center">0.49</td>
<td align="center">1.1</td>
<td align="center">0.42</td>
<td align="center">0.93</td>
<td align="center">23</td>
<td align="center">0.61</td>
<td align="center">44</td>
<td align="center">15</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.8</td>
<td align="center">1.5</td>
<td align="center">0.53</td>
<td align="center">0.32</td>
<td align="center">1.1</td>
<td align="center">0.56</td>
<td align="center">0.93</td>
<td align="center">19</td>
<td align="center">0.4</td>
<td align="center">39</td>
<td align="center">14</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">Vt_ma1</td>
<td align="center">3.4</td>
<td align="center">1.8</td>
<td align="center">0.53</td>
<td align="center">0.38</td>
<td align="center">1.3</td>
<td align="center">0.93</td>
<td align="center">0.94</td>
<td align="center">18</td>
<td align="center">0.45</td>
<td align="center">37</td>
<td align="center">13</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">Vt_ma1</td>
<td align="center">3.1</td>
<td align="center">1.6</td>
<td align="center">0.51</td>
<td align="center">0.34</td>
<td align="center">1.2</td>
<td align="center">0.62</td>
<td align="center">0.93</td>
<td align="center">18</td>
<td align="center">0.43</td>
<td align="center">38</td>
<td align="center">12</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">Vt_ma1</td>
<td align="center">0.93</td>
<td align="center">0.074</td>
<td align="center">0.079</td>
<td align="center">0.016</td>
<td align="center">-0.33</td>
<td align="center">-0.1</td>
<td align="center">0.62</td>
<td align="center">6.1</td>
<td align="center">0.045</td>
<td align="center">14</td>
<td align="center">-1.1</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">Vt_ma1</td>
<td align="center">1</td>
<td align="center">0.035</td>
<td align="center">0.034</td>
<td align="center">0.0074</td>
<td align="center">0.013</td>
<td align="center">-1.2</td>
<td align="center">0.86</td>
<td align="center">1.6</td>
<td align="center">0.013</td>
<td align="center">3.5</td>
<td align="center">0.077</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">Vt_ma1</td>
<td align="center">3</td>
<td align="center">1.6</td>
<td align="center">0.52</td>
<td align="center">0.33</td>
<td align="center">1.1</td>
<td align="center">0.6</td>
<td align="center">0.93</td>
<td align="center">18</td>
<td align="center">0.42</td>
<td align="center">39</td>
<td align="center">9.9</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">Vt_ma1</td>
<td align="center">3.2</td>
<td align="center">1.6</td>
<td align="center">0.52</td>
<td align="center">0.35</td>
<td align="center">1.2</td>
<td align="center">0.93</td>
<td align="center">0.94</td>
<td align="center">17</td>
<td align="center">0.39</td>
<td align="center">34</td>
<td align="center">12</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">Vt_ma1</td>
<td align="center">4.1</td>
<td align="center">2.4</td>
<td align="center">0.58</td>
<td align="center">0.51</td>
<td align="center">1.2</td>
<td align="center">0.57</td>
<td align="center">0.93</td>
<td align="center">21</td>
<td align="center">0.61</td>
<td align="center">41</td>
<td align="center">12</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">Vt_ma1</td>
<td align="center">3.1</td>
<td align="center">1.6</td>
<td align="center">0.51</td>
<td align="center">0.33</td>
<td align="center">1.2</td>
<td align="center">0.79</td>
<td align="center">0.93</td>
<td align="center">18</td>
<td align="center">0.42</td>
<td align="center">38</td>
<td align="center">12</td>
</tr>
</tbody>
</table>

### Separated by Genotype

``` r
HABgroup <- demog$Subjname[demog$Genotype=='HAB']
MABgroup <- demog$Subjname[demog$Genotype=='MAB']
```

``` r
trtout <- trtdata %>%
  filter(Subjname %in% HABgroup) %>%
  gather(Region, Binding, -Subjname, -PETNo, -Measure) %>%
  filter(Subjname != 'mahi') %>%
  spread(PETNo, Binding) %>%
  group_by(Region, Measure) %>%
  do(trt = granviller::trt(.$`1`, .$`2`)$tidy) %>%
  unnest() %>%
  arrange(Measure, Region)

pander::pandoc.table(trtout, digits=2, split.tables=Inf, caption='Test-Retest Analysis for HABs')
```

<table style="width:100%;">
<caption>Test-Retest Analysis for HABs</caption>
<colgroup>
<col width="7%" />
<col width="9%" />
<col width="6%" />
<col width="7%" />
<col width="7%" />
<col width="7%" />
<col width="7%" />
<col width="9%" />
<col width="6%" />
<col width="6%" />
<col width="7%" />
<col width="5%" />
<col width="13%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Region</th>
<th align="center">Measure</th>
<th align="center">mean</th>
<th align="center">sd</th>
<th align="center">cov</th>
<th align="center">se</th>
<th align="center">skew</th>
<th align="center">kurtosis</th>
<th align="center">icc</th>
<th align="center">aapd</th>
<th align="center">sem</th>
<th align="center">md</th>
<th align="center">avgpercchange</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">SUV_4060</td>
<td align="center">1.1</td>
<td align="center">0.27</td>
<td align="center">0.24</td>
<td align="center">0.086</td>
<td align="center">0.2</td>
<td align="center">-1.3</td>
<td align="center">0.8</td>
<td align="center">13</td>
<td align="center">0.12</td>
<td align="center">29</td>
<td align="center">-6</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">SUV_4060</td>
<td align="center">1.1</td>
<td align="center">0.24</td>
<td align="center">0.22</td>
<td align="center">0.075</td>
<td align="center">0.45</td>
<td align="center">-1.2</td>
<td align="center">0.76</td>
<td align="center">13</td>
<td align="center">0.12</td>
<td align="center">30</td>
<td align="center">-3.8</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">SUV_4060</td>
<td align="center">0.94</td>
<td align="center">0.071</td>
<td align="center">0.075</td>
<td align="center">0.022</td>
<td align="center">-1</td>
<td align="center">-0.27</td>
<td align="center">0.85</td>
<td align="center">3.2</td>
<td align="center">0.027</td>
<td align="center">8.1</td>
<td align="center">2.2</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">SUV_4060</td>
<td align="center">1</td>
<td align="center">0.019</td>
<td align="center">0.019</td>
<td align="center">0.0062</td>
<td align="center">0.61</td>
<td align="center">-0.95</td>
<td align="center">0.6</td>
<td align="center">1.3</td>
<td align="center">0.012</td>
<td align="center">3.2</td>
<td align="center">0.26</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">SUV_4060</td>
<td align="center">1</td>
<td align="center">0.24</td>
<td align="center">0.23</td>
<td align="center">0.076</td>
<td align="center">0.17</td>
<td align="center">-0.87</td>
<td align="center">0.74</td>
<td align="center">15</td>
<td align="center">0.12</td>
<td align="center">33</td>
<td align="center">-7</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">SUV_4060</td>
<td align="center">1.1</td>
<td align="center">0.23</td>
<td align="center">0.22</td>
<td align="center">0.073</td>
<td align="center">0.43</td>
<td align="center">-0.94</td>
<td align="center">0.77</td>
<td align="center">13</td>
<td align="center">0.11</td>
<td align="center">28</td>
<td align="center">-6.2</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">SUV_4060</td>
<td align="center">1.4</td>
<td align="center">0.32</td>
<td align="center">0.24</td>
<td align="center">0.1</td>
<td align="center">0.35</td>
<td align="center">-0.86</td>
<td align="center">0.78</td>
<td align="center">14</td>
<td align="center">0.15</td>
<td align="center">31</td>
<td align="center">-6.4</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">SUV_4060</td>
<td align="center">1</td>
<td align="center">0.22</td>
<td align="center">0.21</td>
<td align="center">0.069</td>
<td align="center">0.44</td>
<td align="center">-1</td>
<td align="center">0.78</td>
<td align="center">12</td>
<td align="center">0.1</td>
<td align="center">28</td>
<td align="center">-4</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">SUV_Tot</td>
<td align="center">89</td>
<td align="center">15</td>
<td align="center">0.16</td>
<td align="center">4.6</td>
<td align="center">0.44</td>
<td align="center">-1</td>
<td align="center">0.69</td>
<td align="center">11</td>
<td align="center">8.1</td>
<td align="center">25</td>
<td align="center">-7.1</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">SUV_Tot</td>
<td align="center">87</td>
<td align="center">12</td>
<td align="center">0.14</td>
<td align="center">3.9</td>
<td align="center">0.72</td>
<td align="center">-0.51</td>
<td align="center">0.65</td>
<td align="center">9.9</td>
<td align="center">7.3</td>
<td align="center">23</td>
<td align="center">-5.9</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">SUV_Tot</td>
<td align="center">0.99</td>
<td align="center">0.056</td>
<td align="center">0.057</td>
<td align="center">0.018</td>
<td align="center">-0.47</td>
<td align="center">-1.3</td>
<td align="center">0.96</td>
<td align="center">1.3</td>
<td align="center">0.011</td>
<td align="center">3.1</td>
<td align="center">1.3</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">SUV_Tot</td>
<td align="center">1.1</td>
<td align="center">0.02</td>
<td align="center">0.018</td>
<td align="center">0.0062</td>
<td align="center">0.77</td>
<td align="center">-0.59</td>
<td align="center">0.7</td>
<td align="center">1.2</td>
<td align="center">0.011</td>
<td align="center">2.7</td>
<td align="center">-0.19</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">SUV_Tot</td>
<td align="center">85</td>
<td align="center">13</td>
<td align="center">0.15</td>
<td align="center">4</td>
<td align="center">0.49</td>
<td align="center">-0.65</td>
<td align="center">0.62</td>
<td align="center">11</td>
<td align="center">7.7</td>
<td align="center">25</td>
<td align="center">-6.6</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">SUV_Tot</td>
<td align="center">83</td>
<td align="center">13</td>
<td align="center">0.15</td>
<td align="center">4</td>
<td align="center">0.78</td>
<td align="center">-0.46</td>
<td align="center">0.63</td>
<td align="center">11</td>
<td align="center">7.7</td>
<td align="center">26</td>
<td align="center">-7.2</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">SUV_Tot</td>
<td align="center">100</td>
<td align="center">16</td>
<td align="center">0.16</td>
<td align="center">5</td>
<td align="center">0.7</td>
<td align="center">-0.4</td>
<td align="center">0.71</td>
<td align="center">9.6</td>
<td align="center">8.6</td>
<td align="center">24</td>
<td align="center">-5.8</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">SUV_Tot</td>
<td align="center">79</td>
<td align="center">12</td>
<td align="center">0.15</td>
<td align="center">3.7</td>
<td align="center">0.69</td>
<td align="center">-0.48</td>
<td align="center">0.68</td>
<td align="center">10</td>
<td align="center">6.6</td>
<td align="center">23</td>
<td align="center">-5.7</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">Vt_2tcm</td>
<td align="center">4</td>
<td align="center">1.8</td>
<td align="center">0.44</td>
<td align="center">0.56</td>
<td align="center">0.86</td>
<td align="center">-0.99</td>
<td align="center">0.91</td>
<td align="center">19</td>
<td align="center">0.54</td>
<td align="center">37</td>
<td align="center">7.2</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">Vt_2tcm</td>
<td align="center">3.9</td>
<td align="center">1.6</td>
<td align="center">0.42</td>
<td align="center">0.51</td>
<td align="center">0.73</td>
<td align="center">-1.1</td>
<td align="center">0.89</td>
<td align="center">21</td>
<td align="center">0.54</td>
<td align="center">39</td>
<td align="center">8.9</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">Vt_2tcm</td>
<td align="center">0.97</td>
<td align="center">0.079</td>
<td align="center">0.081</td>
<td align="center">0.025</td>
<td align="center">-0.6</td>
<td align="center">-1.1</td>
<td align="center">0.87</td>
<td align="center">3.4</td>
<td align="center">0.028</td>
<td align="center">7.9</td>
<td align="center">1.8</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">Vt_2tcm</td>
<td align="center">1.1</td>
<td align="center">0.031</td>
<td align="center">0.03</td>
<td align="center">0.0099</td>
<td align="center">-0.24</td>
<td align="center">-1.2</td>
<td align="center">0.33</td>
<td align="center">3.1</td>
<td align="center">0.026</td>
<td align="center">6.8</td>
<td align="center">-1</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">Vt_2tcm</td>
<td align="center">3.6</td>
<td align="center">1.6</td>
<td align="center">0.45</td>
<td align="center">0.52</td>
<td align="center">0.83</td>
<td align="center">-0.98</td>
<td align="center">0.92</td>
<td align="center">19</td>
<td align="center">0.46</td>
<td align="center">35</td>
<td align="center">7.9</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">Vt_2tcm</td>
<td align="center">3.9</td>
<td align="center">1.7</td>
<td align="center">0.44</td>
<td align="center">0.54</td>
<td align="center">0.89</td>
<td align="center">-0.92</td>
<td align="center">0.92</td>
<td align="center">18</td>
<td align="center">0.48</td>
<td align="center">34</td>
<td align="center">7.6</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">Vt_2tcm</td>
<td align="center">5.1</td>
<td align="center">2.5</td>
<td align="center">0.49</td>
<td align="center">0.79</td>
<td align="center">0.75</td>
<td align="center">-1</td>
<td align="center">0.91</td>
<td align="center">24</td>
<td align="center">0.76</td>
<td align="center">41</td>
<td align="center">7.5</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">Vt_2tcm</td>
<td align="center">3.7</td>
<td align="center">1.6</td>
<td align="center">0.44</td>
<td align="center">0.51</td>
<td align="center">0.78</td>
<td align="center">-1</td>
<td align="center">0.91</td>
<td align="center">20</td>
<td align="center">0.5</td>
<td align="center">37</td>
<td align="center">9.9</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">Vt_ma1</td>
<td align="center">4.5</td>
<td align="center">1.9</td>
<td align="center">0.43</td>
<td align="center">0.61</td>
<td align="center">0.92</td>
<td align="center">-0.92</td>
<td align="center">0.92</td>
<td align="center">17</td>
<td align="center">0.55</td>
<td align="center">34</td>
<td align="center">7</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">Vt_ma1</td>
<td align="center">4.1</td>
<td align="center">1.7</td>
<td align="center">0.42</td>
<td align="center">0.55</td>
<td align="center">0.79</td>
<td align="center">-1.1</td>
<td align="center">0.9</td>
<td align="center">20</td>
<td align="center">0.54</td>
<td align="center">37</td>
<td align="center">7.9</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">Vt_ma1</td>
<td align="center">0.92</td>
<td align="center">0.078</td>
<td align="center">0.085</td>
<td align="center">0.025</td>
<td align="center">-0.86</td>
<td align="center">-0.47</td>
<td align="center">0.79</td>
<td align="center">4.8</td>
<td align="center">0.036</td>
<td align="center">11</td>
<td align="center">0.94</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">Vt_ma1</td>
<td align="center">1</td>
<td align="center">0.026</td>
<td align="center">0.026</td>
<td align="center">0.0084</td>
<td align="center">0.53</td>
<td align="center">-1.2</td>
<td align="center">0.79</td>
<td align="center">1.6</td>
<td align="center">0.012</td>
<td align="center">3.3</td>
<td align="center">0.23</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">Vt_ma1</td>
<td align="center">3.9</td>
<td align="center">1.7</td>
<td align="center">0.43</td>
<td align="center">0.53</td>
<td align="center">0.75</td>
<td align="center">-1</td>
<td align="center">0.91</td>
<td align="center">22</td>
<td align="center">0.52</td>
<td align="center">36</td>
<td align="center">3.8</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">Vt_ma1</td>
<td align="center">4.1</td>
<td align="center">1.8</td>
<td align="center">0.43</td>
<td align="center">0.56</td>
<td align="center">0.93</td>
<td align="center">-0.87</td>
<td align="center">0.93</td>
<td align="center">17</td>
<td align="center">0.47</td>
<td align="center">32</td>
<td align="center">5.6</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">Vt_ma1</td>
<td align="center">5.4</td>
<td align="center">2.6</td>
<td align="center">0.48</td>
<td align="center">0.83</td>
<td align="center">0.76</td>
<td align="center">-1</td>
<td align="center">0.91</td>
<td align="center">24</td>
<td align="center">0.77</td>
<td align="center">39</td>
<td align="center">4.1</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">Vt_ma1</td>
<td align="center">4</td>
<td align="center">1.7</td>
<td align="center">0.43</td>
<td align="center">0.54</td>
<td align="center">0.8</td>
<td align="center">-1</td>
<td align="center">0.91</td>
<td align="center">20</td>
<td align="center">0.53</td>
<td align="center">37</td>
<td align="center">7.7</td>
</tr>
</tbody>
</table>

``` r
trtout <- trtdata %>%
  filter(Subjname %in% MABgroup) %>%
  gather(Region, Binding, -Subjname, -PETNo, -Measure) %>%
  filter(Subjname != 'mahi') %>%
  spread(PETNo, Binding) %>%
  group_by(Region, Measure) %>%
  do(trt = granviller::trt(.$`1`, .$`2`)$tidy) %>%
  unnest() %>%
  arrange(Measure, Region)

pander::pandoc.table(trtout, digits=2, split.tables=Inf, caption='Test-Retest Analysis for MABs')
```

<table style="width:100%;">
<caption>Test-Retest Analysis for MABs</caption>
<colgroup>
<col width="7%" />
<col width="9%" />
<col width="6%" />
<col width="7%" />
<col width="7%" />
<col width="7%" />
<col width="7%" />
<col width="9%" />
<col width="6%" />
<col width="6%" />
<col width="7%" />
<col width="5%" />
<col width="13%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Region</th>
<th align="center">Measure</th>
<th align="center">mean</th>
<th align="center">sd</th>
<th align="center">cov</th>
<th align="center">se</th>
<th align="center">skew</th>
<th align="center">kurtosis</th>
<th align="center">icc</th>
<th align="center">aapd</th>
<th align="center">sem</th>
<th align="center">md</th>
<th align="center">avgpercchange</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">SUV_4060</td>
<td align="center">0.84</td>
<td align="center">0.25</td>
<td align="center">0.3</td>
<td align="center">0.072</td>
<td align="center">1</td>
<td align="center">-0.34</td>
<td align="center">0.95</td>
<td align="center">7.3</td>
<td align="center">0.053</td>
<td align="center">17</td>
<td align="center">-1.1</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">SUV_4060</td>
<td align="center">0.8</td>
<td align="center">0.24</td>
<td align="center">0.31</td>
<td align="center">0.071</td>
<td align="center">1.2</td>
<td align="center">0.077</td>
<td align="center">0.91</td>
<td align="center">13</td>
<td align="center">0.073</td>
<td align="center">25</td>
<td align="center">-4.8</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">SUV_4060</td>
<td align="center">0.95</td>
<td align="center">0.05</td>
<td align="center">0.052</td>
<td align="center">0.014</td>
<td align="center">0.38</td>
<td align="center">-1.4</td>
<td align="center">0.32</td>
<td align="center">5.4</td>
<td align="center">0.041</td>
<td align="center">12</td>
<td align="center">-3.7</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">SUV_4060</td>
<td align="center">1</td>
<td align="center">0.041</td>
<td align="center">0.04</td>
<td align="center">0.012</td>
<td align="center">0.59</td>
<td align="center">-1</td>
<td align="center">0.89</td>
<td align="center">1.6</td>
<td align="center">0.013</td>
<td align="center">3.6</td>
<td align="center">0.23</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">SUV_4060</td>
<td align="center">0.77</td>
<td align="center">0.24</td>
<td align="center">0.31</td>
<td align="center">0.069</td>
<td align="center">0.94</td>
<td align="center">-0.3</td>
<td align="center">0.87</td>
<td align="center">16</td>
<td align="center">0.086</td>
<td align="center">31</td>
<td align="center">-4.7</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">SUV_4060</td>
<td align="center">0.8</td>
<td align="center">0.24</td>
<td align="center">0.3</td>
<td align="center">0.069</td>
<td align="center">0.98</td>
<td align="center">-0.32</td>
<td align="center">0.93</td>
<td align="center">9.9</td>
<td align="center">0.062</td>
<td align="center">21</td>
<td align="center">-2.9</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">SUV_4060</td>
<td align="center">1</td>
<td align="center">0.36</td>
<td align="center">0.35</td>
<td align="center">0.1</td>
<td align="center">0.92</td>
<td align="center">-0.32</td>
<td align="center">0.9</td>
<td align="center">16</td>
<td align="center">0.11</td>
<td align="center">30</td>
<td align="center">-3.5</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">SUV_4060</td>
<td align="center">0.79</td>
<td align="center">0.22</td>
<td align="center">0.27</td>
<td align="center">0.062</td>
<td align="center">0.99</td>
<td align="center">-0.12</td>
<td align="center">0.87</td>
<td align="center">13</td>
<td align="center">0.077</td>
<td align="center">27</td>
<td align="center">-5.1</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">SUV_Tot</td>
<td align="center">72</td>
<td align="center">16</td>
<td align="center">0.21</td>
<td align="center">4.5</td>
<td align="center">1.1</td>
<td align="center">-0.021</td>
<td align="center">0.86</td>
<td align="center">11</td>
<td align="center">5.8</td>
<td align="center">22</td>
<td align="center">-4.3</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">SUV_Tot</td>
<td align="center">71</td>
<td align="center">16</td>
<td align="center">0.23</td>
<td align="center">4.8</td>
<td align="center">1.1</td>
<td align="center">0.18</td>
<td align="center">0.8</td>
<td align="center">14</td>
<td align="center">7.4</td>
<td align="center">29</td>
<td align="center">-6.4</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">SUV_Tot</td>
<td align="center">0.98</td>
<td align="center">0.037</td>
<td align="center">0.038</td>
<td align="center">0.011</td>
<td align="center">0.086</td>
<td align="center">-1.6</td>
<td align="center">0.45</td>
<td align="center">3.5</td>
<td align="center">0.027</td>
<td align="center">7.7</td>
<td align="center">-2.1</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">SUV_Tot</td>
<td align="center">1.1</td>
<td align="center">0.038</td>
<td align="center">0.036</td>
<td align="center">0.011</td>
<td align="center">0.39</td>
<td align="center">-1.5</td>
<td align="center">0.89</td>
<td align="center">1.4</td>
<td align="center">0.013</td>
<td align="center">3.4</td>
<td align="center">0.3</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">SUV_Tot</td>
<td align="center">69</td>
<td align="center">17</td>
<td align="center">0.25</td>
<td align="center">4.9</td>
<td align="center">0.74</td>
<td align="center">-0.45</td>
<td align="center">0.77</td>
<td align="center">17</td>
<td align="center">8.1</td>
<td align="center">33</td>
<td align="center">-6</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">SUV_Tot</td>
<td align="center">69</td>
<td align="center">15</td>
<td align="center">0.22</td>
<td align="center">4.3</td>
<td align="center">0.89</td>
<td align="center">-0.23</td>
<td align="center">0.83</td>
<td align="center">12</td>
<td align="center">6.1</td>
<td align="center">25</td>
<td align="center">-5.6</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">SUV_Tot</td>
<td align="center">86</td>
<td align="center">21</td>
<td align="center">0.24</td>
<td align="center">6</td>
<td align="center">0.78</td>
<td align="center">-0.19</td>
<td align="center">0.73</td>
<td align="center">17</td>
<td align="center">11</td>
<td align="center">35</td>
<td align="center">-7</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">SUV_Tot</td>
<td align="center">67</td>
<td align="center">14</td>
<td align="center">0.21</td>
<td align="center">4</td>
<td align="center">0.82</td>
<td align="center">-0.15</td>
<td align="center">0.74</td>
<td align="center">15</td>
<td align="center">7.1</td>
<td align="center">30</td>
<td align="center">-6.7</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.2</td>
<td align="center">1</td>
<td align="center">0.48</td>
<td align="center">0.3</td>
<td align="center">1</td>
<td align="center">-0.24</td>
<td align="center">0.9</td>
<td align="center">23</td>
<td align="center">0.33</td>
<td align="center">42</td>
<td align="center">23</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.2</td>
<td align="center">0.99</td>
<td align="center">0.46</td>
<td align="center">0.29</td>
<td align="center">1</td>
<td align="center">-0.31</td>
<td align="center">0.93</td>
<td align="center">17</td>
<td align="center">0.26</td>
<td align="center">33</td>
<td align="center">17</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">Vt_2tcm</td>
<td align="center">0.99</td>
<td align="center">0.065</td>
<td align="center">0.066</td>
<td align="center">0.019</td>
<td align="center">0.66</td>
<td align="center">-1.4</td>
<td align="center">0.17</td>
<td align="center">5.9</td>
<td align="center">0.06</td>
<td align="center">17</td>
<td align="center">-5.5</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">Vt_2tcm</td>
<td align="center">1</td>
<td align="center">0.038</td>
<td align="center">0.038</td>
<td align="center">0.011</td>
<td align="center">0.19</td>
<td align="center">-1.6</td>
<td align="center">0.56</td>
<td align="center">2.9</td>
<td align="center">0.026</td>
<td align="center">6.9</td>
<td align="center">-0.95</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">Vt_2tcm</td>
<td align="center">2</td>
<td align="center">0.8</td>
<td align="center">0.4</td>
<td align="center">0.23</td>
<td align="center">0.54</td>
<td align="center">-0.76</td>
<td align="center">0.86</td>
<td align="center">18</td>
<td align="center">0.3</td>
<td align="center">41</td>
<td align="center">16</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.1</td>
<td align="center">1</td>
<td align="center">0.47</td>
<td align="center">0.29</td>
<td align="center">0.86</td>
<td align="center">-0.51</td>
<td align="center">0.92</td>
<td align="center">21</td>
<td align="center">0.28</td>
<td align="center">36</td>
<td align="center">21</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.8</td>
<td align="center">1.5</td>
<td align="center">0.54</td>
<td align="center">0.44</td>
<td align="center">1</td>
<td align="center">-0.29</td>
<td align="center">0.93</td>
<td align="center">22</td>
<td align="center">0.41</td>
<td align="center">40</td>
<td align="center">22</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">Vt_2tcm</td>
<td align="center">2.1</td>
<td align="center">0.94</td>
<td align="center">0.45</td>
<td align="center">0.27</td>
<td align="center">0.87</td>
<td align="center">-0.57</td>
<td align="center">0.92</td>
<td align="center">18</td>
<td align="center">0.27</td>
<td align="center">36</td>
<td align="center">18</td>
</tr>
<tr class="odd">
<td align="center">CBL</td>
<td align="center">Vt_ma1</td>
<td align="center">2.5</td>
<td align="center">1.1</td>
<td align="center">0.43</td>
<td align="center">0.31</td>
<td align="center">0.88</td>
<td align="center">-0.37</td>
<td align="center">0.92</td>
<td align="center">19</td>
<td align="center">0.31</td>
<td align="center">34</td>
<td align="center">19</td>
</tr>
<tr class="even">
<td align="center">FC</td>
<td align="center">Vt_ma1</td>
<td align="center">2.3</td>
<td align="center">0.95</td>
<td align="center">0.41</td>
<td align="center">0.28</td>
<td align="center">0.99</td>
<td align="center">-0.33</td>
<td align="center">0.92</td>
<td align="center">16</td>
<td align="center">0.28</td>
<td align="center">33</td>
<td align="center">16</td>
</tr>
<tr class="odd">
<td align="center">FC_CBL</td>
<td align="center">Vt_ma1</td>
<td align="center">0.94</td>
<td align="center">0.071</td>
<td align="center">0.075</td>
<td align="center">0.021</td>
<td align="center">0.37</td>
<td align="center">-1.4</td>
<td align="center">0.48</td>
<td align="center">7.1</td>
<td align="center">0.051</td>
<td align="center">15</td>
<td align="center">-2.7</td>
</tr>
<tr class="even">
<td align="center">FC_WB</td>
<td align="center">Vt_ma1</td>
<td align="center">1</td>
<td align="center">0.04</td>
<td align="center">0.039</td>
<td align="center">0.012</td>
<td align="center">0.19</td>
<td align="center">-1.7</td>
<td align="center">0.89</td>
<td align="center">1.6</td>
<td align="center">0.013</td>
<td align="center">3.6</td>
<td align="center">-0.049</td>
</tr>
<tr class="odd">
<td align="center">STR</td>
<td align="center">Vt_ma1</td>
<td align="center">2.3</td>
<td align="center">0.95</td>
<td align="center">0.42</td>
<td align="center">0.27</td>
<td align="center">0.91</td>
<td align="center">-0.31</td>
<td align="center">0.9</td>
<td align="center">16</td>
<td align="center">0.3</td>
<td align="center">36</td>
<td align="center">15</td>
</tr>
<tr class="even">
<td align="center">TC</td>
<td align="center">Vt_ma1</td>
<td align="center">2.3</td>
<td align="center">0.99</td>
<td align="center">0.42</td>
<td align="center">0.29</td>
<td align="center">0.83</td>
<td align="center">-0.45</td>
<td align="center">0.93</td>
<td align="center">16</td>
<td align="center">0.27</td>
<td align="center">32</td>
<td align="center">16</td>
</tr>
<tr class="odd">
<td align="center">THA</td>
<td align="center">Vt_ma1</td>
<td align="center">3</td>
<td align="center">1.5</td>
<td align="center">0.51</td>
<td align="center">0.43</td>
<td align="center">0.99</td>
<td align="center">-0.37</td>
<td align="center">0.93</td>
<td align="center">19</td>
<td align="center">0.38</td>
<td align="center">36</td>
<td align="center">19</td>
</tr>
<tr class="even">
<td align="center">WB</td>
<td align="center">Vt_ma1</td>
<td align="center">2.3</td>
<td align="center">0.91</td>
<td align="center">0.4</td>
<td align="center">0.26</td>
<td align="center">0.82</td>
<td align="center">-0.59</td>
<td align="center">0.91</td>
<td align="center">16</td>
<td align="center">0.27</td>
<td align="center">33</td>
<td align="center">16</td>
</tr>
</tbody>
</table>

Prepare data for Plots
----------------------

``` r
plotdat <- data.frame(PET = datdf$PET, Subjname = datdf$Subjname,
                     Genotype=datdf$Genotype, PETNo=datdf$PETNo,
                     SUV_FC_WB = SUV_tot$FC_WB,
                     SUV_FC_CBL = SUV_tot$FC_CBL, 
                     SUV4060_FC_WB = SUV_4060$FC_WB,
                     SUV4060_FC_CBL = SUV_4060$FC_CBL,
                     DVR_FC_WB = VT_2tcm$FC_WB,
                     DVR_FC_CBL = VT_2tcm$FC_CBL,
                     VT_FC = VT_2tcm$FC,
                     VT_FC_MA1 = VT_MA1$FC,
                     SUV_FC = SUV_4060$FC) %>%
  filter(Subjname != 'mahi')
```

Calculate the correlations between metrics
------------------------------------------

``` r
corout <- plotdat %>%
  gather(Measure, Binding, -PET, -Genotype, -PETNo, -VT_FC, -Subjname) %>%
  group_by(Measure, Genotype) %>%
  dplyr::summarise(r2=cor(Binding, VT_FC)^2) %>%
  arrange(Genotype, r2)

pander::pandoc.table(corout, digits=2, split.tables=Inf)
```

<table style="width:54%;">
<colgroup>
<col width="23%" />
<col width="15%" />
<col width="15%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Measure</th>
<th align="center">Genotype</th>
<th align="center">r2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">SUV_FC_CBL</td>
<td align="center">HAB</td>
<td align="center">0.00096</td>
</tr>
<tr class="even">
<td align="center">SUV4060_FC_CBL</td>
<td align="center">HAB</td>
<td align="center">0.011</td>
</tr>
<tr class="odd">
<td align="center">DVR_FC_CBL</td>
<td align="center">HAB</td>
<td align="center">0.016</td>
</tr>
<tr class="even">
<td align="center">SUV4060_FC_WB</td>
<td align="center">HAB</td>
<td align="center">0.024</td>
</tr>
<tr class="odd">
<td align="center">SUV_FC_WB</td>
<td align="center">HAB</td>
<td align="center">0.15</td>
</tr>
<tr class="even">
<td align="center">DVR_FC_WB</td>
<td align="center">HAB</td>
<td align="center">0.31</td>
</tr>
<tr class="odd">
<td align="center">SUV_FC</td>
<td align="center">HAB</td>
<td align="center">0.69</td>
</tr>
<tr class="even">
<td align="center">VT_FC_MA1</td>
<td align="center">HAB</td>
<td align="center">0.99</td>
</tr>
<tr class="odd">
<td align="center">SUV4060_FC_CBL</td>
<td align="center">MAB</td>
<td align="center">0.0034</td>
</tr>
<tr class="even">
<td align="center">DVR_FC_WB</td>
<td align="center">MAB</td>
<td align="center">0.0055</td>
</tr>
<tr class="odd">
<td align="center">SUV_FC_CBL</td>
<td align="center">MAB</td>
<td align="center">0.028</td>
</tr>
<tr class="even">
<td align="center">DVR_FC_CBL</td>
<td align="center">MAB</td>
<td align="center">0.034</td>
</tr>
<tr class="odd">
<td align="center">SUV_FC_WB</td>
<td align="center">MAB</td>
<td align="center">0.21</td>
</tr>
<tr class="even">
<td align="center">SUV4060_FC_WB</td>
<td align="center">MAB</td>
<td align="center">0.33</td>
</tr>
<tr class="odd">
<td align="center">SUV_FC</td>
<td align="center">MAB</td>
<td align="center">0.86</td>
</tr>
<tr class="even">
<td align="center">VT_FC_MA1</td>
<td align="center">MAB</td>
<td align="center">0.99</td>
</tr>
</tbody>
</table>

Make interregional correlation plots for V<sub>T</sub>
------------------------------------------------------

``` r
library(corrplot)

col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))

jpeg(filename = 'CorrelationPlots.jpg', width = 9, height=8, res = 800, units = 'in')
par(mfrow=c(2,2))
VT_2tcm %>%
  filter(Subjname %in% HABgroup) %>%
  select(CBL:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 col=col2(200), diag='n',
                 number.digits = 2, title=expression(HAB ~ V[T] ~ Correlations),
                 mar=c(0,0,1,0))

VT_2tcm %>%
  filter(Subjname %in% MABgroup) %>%
  select(CBL:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 col=col2(200), diag='n',
                 number.digits = 2, , title=expression(MAB ~ V[T] ~ Correlations),
                 mar=c(0,0,1,0))

SUV_4060 %>%
  filter(Subjname %in% HABgroup) %>%
  select(CBL:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 col=col2(200), diag='n',
                 number.digits = 2, title=expression(HAB ~ SUV ~ Correlations),
                 mar=c(0,0,1,0))

SUV_4060 %>%
  filter(Subjname %in% MABgroup) %>%
  select(CBL:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 col=col2(200), diag='n',
                 number.digits = 2, title=expression(MAB ~ SUV ~ Correlations),
                 mar=c(0,0,1,0))
dev.off()
```

    ## png 
    ##   2

![Interregional correlations](CorrelationPlots.jpg)

Perform a PCA
-------------

Note that I perform scaling both within genotype and region. Thus HABs and MABs can be combined in the same PCA analysis.

``` r
pcadat <- longdat %>%
  ungroup %>%
  inner_join(select(datdf, PET, Genotype)) %>%
  select(PET, Subjname, PETNo, Genotype, Region, Vt_2tcm) %>%
  group_by(Region, Genotype) %>%
  mutate(Vt_2tcm.z = scale(Vt_2tcm)) %>%
  ungroup() %>%
  select(PET, Subjname, PETNo, Genotype, Region, Vt_2tcm.z) %>%
  filter(PET != 'mahi_2')
```

    ## Joining, by = "PET"

### On all regions

``` r
pca1 <- pcadat %>%
  filter(PETNo==1) %>%
  spread(Region, Vt_2tcm.z) %>%
  select(CBL:WB) %>%
  prcomp() %>%
  summary()

pander::pandoc.table(pca1$importance, digits=3, caption = "PCA for PET1: All Regions")
```

<table style="width:100%;">
<caption>PCA for PET1: All Regions</caption>
<colgroup>
<col width="33%" />
<col width="9%" />
<col width="11%" />
<col width="10%" />
<col width="11%" />
<col width="11%" />
<col width="11%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"> </th>
<th align="center">PC1</th>
<th align="center">PC2</th>
<th align="center">PC3</th>
<th align="center">PC4</th>
<th align="center">PC5</th>
<th align="center">PC6</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong>Standard deviation</strong></td>
<td align="center">2.37</td>
<td align="center">0.227</td>
<td align="center">0.122</td>
<td align="center">0.0832</td>
<td align="center">0.0294</td>
<td align="center">0.0169</td>
</tr>
<tr class="even">
<td align="center"><strong>Proportion of Variance</strong></td>
<td align="center">0.987</td>
<td align="center">0.00902</td>
<td align="center">0.0026</td>
<td align="center">0.00122</td>
<td align="center">0.00015</td>
<td align="center">5e-05</td>
</tr>
<tr class="odd">
<td align="center"><strong>Cumulative Proportion</strong></td>
<td align="center">0.987</td>
<td align="center">0.996</td>
<td align="center">0.999</td>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">1</td>
</tr>
</tbody>
</table>

``` r
pca2 <- pcadat %>%
  filter(PETNo==2) %>%
  spread(Region, Vt_2tcm.z) %>%
  select(CBL:WB) %>%
  prcomp() %>%
  summary()

pander::pandoc.table(pca2$importance, digits=3, caption = "PCA for PET2: All Regions", split.tables=Inf)
```

<table style="width:100%;">
<caption>PCA for PET2: All Regions</caption>
<colgroup>
<col width="33%" />
<col width="9%" />
<col width="11%" />
<col width="11%" />
<col width="11%" />
<col width="11%" />
<col width="11%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"> </th>
<th align="center">PC1</th>
<th align="center">PC2</th>
<th align="center">PC3</th>
<th align="center">PC4</th>
<th align="center">PC5</th>
<th align="center">PC6</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong>Standard deviation</strong></td>
<td align="center">2.49</td>
<td align="center">0.144</td>
<td align="center">0.0851</td>
<td align="center">0.0762</td>
<td align="center">0.0408</td>
<td align="center">0.0281</td>
</tr>
<tr class="even">
<td align="center"><strong>Proportion of Variance</strong></td>
<td align="center">0.994</td>
<td align="center">0.00332</td>
<td align="center">0.00116</td>
<td align="center">0.00093</td>
<td align="center">0.00027</td>
<td align="center">0.00013</td>
</tr>
<tr class="odd">
<td align="center"><strong>Cumulative Proportion</strong></td>
<td align="center">0.994</td>
<td align="center">0.998</td>
<td align="center">0.999</td>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">1</td>
</tr>
</tbody>
</table>

### Only FC, WB and CBL

``` r
pca1 <- pcadat %>%
  filter(PETNo==1) %>%
  spread(Region, Vt_2tcm.z) %>%
  select(FC, CBL, WB) %>%
  prcomp() %>%
  summary()

pander::pandoc.table(pca1$importance, digits=3, caption = "PCA for PET1: FC, WB and CBL")
```

<table style="width:79%;">
<caption>PCA for PET1: FC, WB and CBL</caption>
<colgroup>
<col width="40%" />
<col width="11%" />
<col width="13%" />
<col width="13%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"> </th>
<th align="center">PC1</th>
<th align="center">PC2</th>
<th align="center">PC3</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong>Standard deviation</strong></td>
<td align="center">1.68</td>
<td align="center">0.119</td>
<td align="center">0.0249</td>
</tr>
<tr class="even">
<td align="center"><strong>Proportion of Variance</strong></td>
<td align="center">0.995</td>
<td align="center">0.00499</td>
<td align="center">0.00022</td>
</tr>
<tr class="odd">
<td align="center"><strong>Cumulative Proportion</strong></td>
<td align="center">0.995</td>
<td align="center">1</td>
<td align="center">1</td>
</tr>
</tbody>
</table>

``` r
pca2 <- pcadat %>%
  filter(PETNo==2) %>%
  spread(Region, Vt_2tcm.z) %>%
  select(FC, CBL, WB) %>%
  prcomp() %>%
  summary()

pander::pandoc.table(pca2$importance, digits=3, caption = "PCA for PET2: FC, WB and CBL")
```

<table style="width:79%;">
<caption>PCA for PET2: FC, WB and CBL</caption>
<colgroup>
<col width="40%" />
<col width="11%" />
<col width="13%" />
<col width="13%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"> </th>
<th align="center">PC1</th>
<th align="center">PC2</th>
<th align="center">PC3</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong>Standard deviation</strong></td>
<td align="center">1.76</td>
<td align="center">0.0807</td>
<td align="center">0.0534</td>
</tr>
<tr class="even">
<td align="center"><strong>Proportion of Variance</strong></td>
<td align="center">0.997</td>
<td align="center">0.00209</td>
<td align="center">0.00091</td>
</tr>
<tr class="odd">
<td align="center"><strong>Cumulative Proportion</strong></td>
<td align="center">0.997</td>
<td align="center">0.999</td>
<td align="center">1</td>
</tr>
</tbody>
</table>

Make the plot
-------------

``` r
# devtools::install_github("mvuorre/vmisc")
library(vmisc)
library(gridExtra)
library(grid)
library(RColorBrewer)

a <- ggplot(plotdat, aes(x=VT_FC, y=SUV4060_FC_WB, colour=Genotype)) + 
  geom_point() + expand_limits(x=0,y=c(0.95,1.1)) + theme_blog() + 
  labs(x=expression(V[T]),
       y=expression(SUVR[WB])) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.background = element_rect(fill = "white")) +
  scale_color_brewer(palette = 'Set1') + geom_smooth(method="lm", se=F)

b <- ggplot(plotdat, aes(x=VT_FC, y=SUV4060_FC_CBL, colour=Genotype)) + 
  geom_point() + expand_limits(x=0,y=c(0.8,1.1)) + theme_blog() + 
  labs(x=expression(V[T]),
       y=expression(SUVR[CBL])) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks=seq(0.8,1.1,by=0.1)) +
  scale_color_brewer(palette = 'Set1') + geom_smooth(method="lm", se=F)

c <- ggplot(plotdat, aes(x=VT_FC, y=DVR_FC_WB, colour=Genotype)) + 
  geom_point() + expand_limits(x=0,y=c(0.95,1.1)) + theme_blog() + 
  labs(x=expression(paste('Frontal Cortex ',V[T])),
       y=expression(DVR[WB])) +
  theme(axis.title.y=element_blank(),
        plot.background = element_rect(fill = "white")) +
  scale_color_brewer(palette = 'Set1') + geom_smooth(method="lm", se=F)

d <- ggplot(plotdat, aes(x=VT_FC, y=DVR_FC_CBL, colour=Genotype)) + 
  geom_point() + expand_limits(x=0,y=c(0.8,1.1)) + theme_blog() + 
  labs(x=expression(paste('Frontal Cortex ',V[T])),
       y=expression(DVR[CBL])) +
  theme(axis.title.y=element_blank(),
        plot.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks=seq(0.8,1.1,by=0.1)) +
  scale_color_brewer(palette = 'Set1') + geom_smooth(method="lm", se=F)

e = textGrob(label = 'SUVR')
f = textGrob(label = 'DVR')
g = textGrob(label = '     Whole Brain')
h = textGrob(label = '     Cerebellum')

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

a2 <- ggplot(plotdat, aes(x=VT_FC, y=SUV4060_FC_WB, colour=Genotype)) + 
  geom_point() + theme_bw() +
  scale_color_brewer('Genotype', palette = 'Set1')

genlegend <- g_legend(a2)

layout <- rbind(c(NA, 7, 8, 9 ),
                c(1, 2, 3, 9),
                c(4, 5, 6, 9))

gridExtra::grid.arrange(grobs=list(e,a,b,f,c,d,g,h, genlegend), layout_matrix=layout, widths=c(1,4,4,2), heights=c(1,4,5))
```

![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
# svglite::svglite('abstractFig.svg', width = 7.8, height = 6)
# gridExtra::grid.arrange(grobs=list(e,a,b,f,c,d,g,h, genlegend), layout_matrix=layout)
# dev.off()
```

``` r
alabel <- textGrob(label = expression(paste('MA1 FC ',V[T])))
blabel <- textGrob(label = expression(paste('FC ',SUV[40-60])))

Vta <- ggplot(plotdat, aes(x=VT_FC, y=VT_FC_MA1, colour=Genotype)) + 
  geom_point() + theme_blog() + 
  labs(x=expression(paste('2TCM Frontal Cortex ',V[T])),
       y=expression(paste('MA1 Frontal Cortex ',V[T]))) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks=seq(1,8,by=0.5)) +
  scale_color_brewer(palette = 'Set1') + geom_smooth(method="lm", se=F) +
  geom_abline(slope = 1, linetype='dashed')

Vtb <- ggplot(plotdat, aes(x=VT_FC, y=SUV_FC, colour=Genotype)) + 
  geom_point() + theme_blog() + 
  labs(x=expression(paste('2TCM Frontal Cortex ',V[T])),
       y=expression(paste('Frontal Cortex ',SUVV[40-60]))) +
  theme(axis.title.y=element_blank(),
        plot.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks=seq(0.5,1.6,by=0.1)) +
  scale_color_brewer(palette = 'Set1') + geom_smooth(method="lm", se=F)

layout2 <- rbind(c(3, 1, 5),
                 c(4, 2, 5))

grid.arrange(grobs=list(Vta, Vtb, alabel, blabel, genlegend), layout_matrix=layout2, widths=c(2,4,2), heights=c(4,5))
```

![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
longdat_showFits <- longdat %>%
  filter(Region %in% c('FC', 'WB', 'CBL'))

delayFits = map(longdat_showFits$WB_delay[longdat_showFits$Region=='WB'], ~plot_inptac_fit(.x) + ggtitle('Delay'))
allFits_2tcm = map2(longdat_showFits$fit_2tcm, longdat_showFits$Region, ~plot_kinfit(.x, roiname=.y))

PETs <- unique(longdat_showFits$PET)

delayFits <- data.frame(PET = PETs) %>%
  mutate(Fit = delayFits,
         Plot = 'Delay')

allFits <- data.frame(PET = rep(PETs, each=3)) %>%
  mutate(Fit = allFits_2tcm,
         Plot = 'Fit',
         PET = as.character(PET)) %>%
  bind_rows(delayFits) %>%
  arrange(PET, Plot)

allFits_excluded <- allFits %>%
  filter(grepl(PET, pattern='mahi'))

allFits <- allFits %>%
  filter(!grepl(PET, pattern='mahi'))

fitLabels <- unique(allFits$PET)

marrangeGrob(allFits$Fit, nrow=2, ncol=2, top=quote(paste('PET: ', PETs[g])))
```

![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-1.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-2.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-3.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-4.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-5.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-6.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-7.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-8.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-9.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-10.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-11.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-12.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-13.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-14.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-15.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-16.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-17.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-18.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-19.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-20.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-21.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-25-22.png)

And the excluded participant:

``` r
fitLabels_excl <- unique(allFits_excluded$PET)
marrangeGrob(allFits_excluded$Fit, nrow=2, ncol=2, top=quote(paste('PET: ', fitLabels_excl[g])))
```

![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-26-1.png)![](ModellingAnalysis_files/figure-markdown_github/unnamed-chunk-26-2.png)
