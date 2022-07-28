

# Cleaning my repositories of Dolo.jl

As of 14.07, the folder GitHub is my "large" working directory for github, all the repository I would clone would end up in this folder, and if I fork a repository, I would put it here.

Now GitHub has 2 folders of interest

- Dolo.jl, which is the cloned repository, it has the latest updates, however it has been wrongly linked to a fork of Dolo.py (which I do not remember forking!)

- Pablo winant internship includes a version of Dolo.jl (that I intentionally forked! and updated) but this version is not correct since it does not contain a file .git 

- the way I updated Dolo.jl: I did a pull changes from Econforge/Dolo.jl/master to eliottmlg/Dolo.jl/master

## what i want to do 

- either delete the dolo.py repo on my github, hoping it will cut the link with Dolo.jl cloned (while not deleting the local folder either!) and then create a new repository from a local file 
- but this would the wrong procedure since it would not be an official fork
- so maybe we cloned it just bc it was quicker to check that indeed the functions work 

- second option is to fork the current version of Econforge/Dolo.jl/master, this would be an official fork, and then I can add the folder eliottexp in it, and pull requests.
- this would require that I first save folder eliottexp in julia_training--this is done already
- also requires that I delete the current fork, so would delete it on github and in GitHub directory.
- need to check everything in Dolo.jl forked (time iteration) to make sure I did not make important changes 
- then delete Dolo.py


- Third option: update Dolo.jl forked
