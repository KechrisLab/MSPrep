
# Updated version of functions from Sean Jacobson

<<<<<<< HEAD

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


=======
Email from 6/29:

> I've been working with MSPrep for a couple of years now, and the original
> version (the version that you have) has quite a few issues with it, to say the
> least. I've done a great deal of work fixing bugs and inconsistincies in the
> code. In order to avoid you repeating a lot of the work I've done, we should
> have some discussions about what I've done and where you can go from there.
> 
> To start, I've attached the file with all of the functions that I've edited.
> Almost every function in the original MSPrep package had major issues and
> needed to be fixed. The two functions I've never changed are graphimputations
> and diagnosticgen, but that is because I never use them. I would be very
> surprised if they didn't have major problems as well.
> 
> In the file I've attached, normdata.new\_kk is the version of the "normdata"
> code in the original package that Dr. Kechris fixed. Her fix allowed the code
> to run without errors, but it didn't fix the code so that it could take
> multiple phenotypes as input. Irritatingly, if you try to input multiple
> phenotypes, the function won't give you an error, but rather use only the first
> phenotype in the input and throw the rest out.
> 
> I created a function normdata.new\_multpheno which can take as input (and
> actually use) (A) Multiple phenotypes (B) No phenotypes and (C) exactly one
> phenotype. The previous function could only do (C). The problem with my
> function is that it only does crmn and ComBat - the other function did more
> than that. I never fixed it for the other types of normalization because I
> wasn't using those. Also, the previous function log-normalized the data in the
> middle of doing the ComBat normalization - I've made that optional, in case the
> user wants to implement a different type of normalization.


Overview:
- resolve differences in normdata functions. Versions: 
  1. `normdata.new_kk` (vers. fixed by Katerina)
  1. `normdata.new_multpheno` (vers. updated by Sean)
  - Also a version in 5/5 version (main working version) and possibly
    Harrison's 7/5 version
>>>>>>> seans_v1.1
