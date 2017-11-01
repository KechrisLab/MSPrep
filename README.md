
MSPrep R package 


Summary from meeting w/ Domink:
- Start with read\_data
- Sean's version (in sean branch, \_sj()) has more code than the one in
  develop/master
- Separate out tasks is the first key thing to do

Notes from meeting w/ Dominik:

  - Dominik uses:
    - read_data_sj()
      - fn doesn't allow more than 3 technical replicates -- should allow Inf
      - fn summarizes, but this should be a separate fn -- probabaly summary()
      - option for whether data is already log-xfrmed or not - 
      - option to load/convert data existing in workign environment
      - ID variable ("lcms...") is hardcoded, shoudln't be
  - batch correction -- 
    - hasn't used
    - should this be done before summaries, imputation, normalization
  - cvmax - 
  - missing data typically represented by 0 or 1 (log(1)=0)


