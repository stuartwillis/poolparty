# PoolParty

### A BASH pipeline to align and analyze paired-end NGS data on a genome assembly.  

Citation: "Micheletti SJ and SR Narum. 2018. Utility of pooled sequencing for association mapping in nonmodel organisms. Molecular Ecology Resources 10.1111/1755-0998.12784"

## Getting Started

 PoolParty is designed to be run on Linux servers. As such, memory and storage may be limiting factors for some systems depending on genome sizes, number of pools, number of SNPs, etc.

 It is highly recommended to run the example files provided in the example directories before analyzing real datasets.  

 PDF files in the example folders contains additional detailed information on the pipeline.  
 
 We stronly recommend the use of linux 'screens' when running the PoolParty modules, as this avoids connection-interription crashes.

 Instructions are provided below for installing dependencies with conda, i.e. without sudo privileges. 
 
 
## Dependencies

PoolParty is designed to be run on Linux (GNU) operating systems. Because it coordinates the execution of multiple packages there are number of dependencies that must be installed prior to running. With the use of diverse packages, the latest versions of Java, Perl, and R must be installed. The required packages for PoolParty are:

### Required package with version at inception 
- Burrows-Wheeler Aligner (BWA; 07.12) - http://bio-bwa.sourceforge.net/  
- Fastqc (0.11.7 ) - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
- samblaster (0.1.24) - https://github.com/GregoryFaust/samblaster  
- samtools (1.5) - http://www.htslib.org/download/  
- bcftools (1.5) - http://www.htslib.org/download/  
- Picard Tools (2.17.11) - http://broadinstitute.github.io/picard/  
- Popoolation2 (1.201) - https://sourceforge.net/p/popoolation2/wiki/Main/  
- BBMap (37.93) - https://sourceforge.net/projects/bbmap/  


### R packages used
If not already installed, PoolParty will attempt to automatically install required R packages. It is recommended to manually install packages beforehand:  

-PPalign: matrixStats, tidyr, stringr, data.table  
-PPstats: reshape, fBasics, ggplot2, RColorBrewer  
-PPanalyze: matrixStats, plyr, stringr, data.table, fBasics, ape, metap

## Installing Conda

#### The basis for this installation is **conda** (Anaconda/Miniconda), which was designed for python but has been extended for a multitude of common programs. Conda is an excellent software manager that provides for software installation without high-level permissions (e.g. root or sudo), or in most cases, worrying about system idiosyncrasies (RHEL vs. Ubuntu, library presence/location, etc.). It is an *"It just works"* provider

There are generally two options for installing with conda; each user has a global conda 'environment' that provides programs in THEIR path (but no one else's); this may still suffer some version conflicts; the other option is to create separate environments where different copies/versions of programs are installed; this prevents version conflicts, but each 'environment' must be activated each time to make those programs accessible in the path (which is simple)

Though most users will not notice a difference, there are also two flavors of conda: Anaconda and Miniconda. The difference is that Anaconda is a much larger installation with many included programs, most of which most users will never use. Miniconda is the light version, installing much faster; we recommend that one below, though either should serve for our purposes.

We install **conda** as follows, beginning in our home directory

> cd ~

> wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

> bash Miniconda3-latest-Linux-x86_64.sh

This will install the latest miniconda software in your home directory. Answer 'yes' as needed. If any problems are encountered, consult the conda help literature. You may need to start a new ssh session for conda to activate.

### Adding conda channels, with priority

Conda works by drawing 'recipes' for installation from 'channels'. These opt-in channels are stored in the `.condarc` file (this hidden file resides in your home directory; find it with `ls -a`). Each time a conda channel is added, it takes priority over the previous ones, so that if an identically-named recipe is found in multiple, it will be pulled from the higher priority channel. We add channels as follows:

> conda config --add channels defaults

> conda config --add channels conda-forge

> conda config --add channels bioconda

> conda config --add channels r

FYI, the 'biobuilds' channel is necessary for Mac users running PoolParty in terminal to install a Mac-compatible samblaster. If that's not you, skip it.

> conda config --add channels biobuilds 

Take a look at your channels

> cat .condarc

to check that channels were added properly. If one wished to change the priority of a channel, just add it again to the top as before, though one can also specify which channel to pull a recipe from a specific channel.

Notice also that the channel priority defaults to "strict". Although we prefer that conda get a package from channels in the order we've specified (in .condarc), we usually prefer that conda solve version conflicts itself rather than failing and reporting why, which it does by choosing from the other channels we have listed. In other words, we want it to be "flexible" in terms of choosing (priority first, then solve conflicts). We can tell it to do that with the following command:

> conda config --set channel_priority flexible

Now take a look at your channels files again, and would should see that "channel priority" is set to "flexible".

> cat .condarc

## Installing The Pipeline and its Dependencies

All the following installation could theoretically be done in the user's global ("base") conda environment. But as described, it's safer to create separate environments for different program dependency suites. Let's create one for PoolParty, and in doing so, install the R environment from the conda-forge environment.

> conda create -n poolparty_env -c conda-forge r-base=4

We have to specify to make sure R 4.x is installed as many packages are now only available for basic install with 4+. And we have to install r-gert for the "devtools" package, which we need to install an older version of metap, because the newer version uses 'qqconf', which Installation may take a while. Once it's finishes, we will have a conda environment named 'poolparty_env' with R installed. This environment will not automatically be active. To activate this environment, we use

> conda activate poolparty_env

or

> source activate poolparty_env

depending on your linux distribution. Once it's active, your prompt should be preceded with something like

`(poolparty_env) [user]$`

To return to your base environment, 'deactivate', 'conda activate', or just log out (but wait until you've installed the remaining dependencies and tested poolparty!) Let's check that R was installed properly:

> which R

should return something like

`~/miniconda3/envs/poolparty_env/bin/R`

If not, check for errors and try again. Now start R

> R

And enter the following installation commands

```
install.packages("BiocManager")
BiocManager::install(c("BiocGenerics","Biobase"))
install.packages(c("survival","MASS"))
BiocManager::install("multtest")
library("multtest")
install.packages(c("metap","ape","matrixStats","fBasics","bibtex","gbRd","Rdpack"))
install.packages(c("ggplot2","RColorBrewer","data.table","tidyr"))
BiocManager::install("qvalue")
quit(save="no")
```

You will need to resolve any installation errors before proceeding. One occasional bug is that the most recent version of 'metap' can't be installed remotely because it now depends on 'qqconf' which also can't be installed because it can't find some FFTW libraries (don't get me started). Alternatively, we can download (outside of R) and older version and install locally, as follows:

> wget https://cran.r-project.org/src/contrib/Archive/metap/metap_1.7.tar.gz

> R

```
install.packages("metap_1.7.tar.gz",repos=NULL, type="source")
library(metap)
quit(save="no")
```

Another occasional bug is that 'multtest', from Bioconductor, can't be installed remotely. Instead, we can also download it outside of R and repeat the installation. Check the multest website (https://www.bioconductor.org/packages/release/bioc/html/multtest.html) to confirm which version is available. 

> wget https://www.bioconductor.org/packages/release/bioc/src/contrib/multtest_2.50.0.tar.gz

> R

```
install.packages("multtest_2.50.0.tar.gz",repos=NULL, type="source")
library("multtest")
quit(save="no")
```

Now let's install the other PoolParty dependencies. These should install from the bioconda channel, which should take precedence in the order of channels (by virtue of being added later), but flexibly added from conda-forge to resolve conflicts.

> conda install bwa fastqc samblaster samtools bcftools picard bbmap perl-app-cpanminus parallel dos2unix

Answer 'yes' as necessary. Now we need to install the perl module *PPanalyze* uses for Fisher's exact test, using CPAN (hence why we installed perl-app-cpanminus)

> cpan Text::NSP::Measures::2D::Fisher::twotailed

Answer yes if asked about auto-configuration. Now make sure bbmap can find its adapter file (confirm they're present, then copy them)

> ls ~/miniconda3/envs/poolparty_env/opt/bbmap*/resources

If this list files rather than throws an error, you're good to copy. If not, find where the adapter files are and change the 'source' path below.

> cp -r ~/miniconda3/envs/poolparty_env/opt/bbmap*/resources ~/miniconda3/envs/poolparty_env/bin

Check that samtools works.

> samtools

It may complain about libraries, which sometimes occurs because R 4+ installs a different version of openssl than samtools expects (1.0.0). If a library can't be found, try linking a similar library (e.g. libcrypto.so, libcrypto.so.3, libcrypto.so.1.1, etc.; AFAIK they all work) in it's place. Modify the names as appropriate in the second command below. Otherwise, if samtools works (lists its options with the above call), ignore this step.

> ls ~/miniconda3/envs/poolparty_env/lib/libcrypto.so*

> ln -s ~/miniconda3/envs/poolparty_env/lib/libcrypto.so ~/miniconda3/envs/poolparty_env/lib/libcrypto.so.1.0.0

> samtools

Finally, Popoolation2 and PoolParty don't actually require any installation (they are just perl and bash scripts), so we can just put them in a special user directory for programs, ~/bin

> cd ~

> mkdir bin

> cd bin

> wget --no-check-certificate https://downloads.sourceforge.net/project/popoolation2/popoolation2_1201.zip

> unzip popoolation2_1201.zip

> rm popoolation2_1201.zip

> wget https://github.com/stuartwillis/poolparty/archive/refs/heads/main.zip

> unzip main.zip

> rm main.zip

> mv poolparty-main poolparty

These now reside at `~/bin`. You shouold specify `~/bin/popoolation2_1201` in the config files for PoolParty as needed. PRO tip: you can add ~/bin to your user $PATH, so that everything in that folder is globally accessible in the $PATH, by adding it to the line in ~/.bash_profile as '$HOME/bin'. 

Now, we will need to specify the path to the jar file for Picard, which we installed with conda. Let's check where it is:

> ls ~/miniconda3/envs/poolparty_env/share/picard*/*.jar

Depending on what this returns, specify `~/miniconda3/envs/poolparty_env/share/picard-2.20.2-0/picard.jar` *or similar* as needed in the config files.

You should now be ready to run PoolParty. Try starting with the examples ;-)

## Troubleshooting

PoolParty is constantly evolving so users may encounter bugs. However, there are common issues that can be avoided  :

1) Permissions: Proper permissions are not only needed to run the PoolParty modules, but also all dependencies. Ensure that your user account has permissions to execute programs and write to the specified output directories.
2) Memory: With increased data comes increased memory usage. If java programs encounter a memory error they will usually spit out an interpretable error. Tune the java memory parameter accordingly.
3) Storage: Large temporary files can fill up smaller hard drives fast. Storage issues generally will have to be resolved with hardware. 
4) Compatibility: PoolParty is POSIX compliant but incompatibilities with specific Linux distributions can still be encountered. Specific formatting of output drives can cause issues, especially if piping is not supported on these drives (mkfifo). Errors associating with drives may require reformatting or diverting output files to a different drive. 

If an issue does not fall within this category, post the error message and explanation to the PoolParty GitHub page. Also, don't forget to check out the example file for more details. 

## Stopping PoolParty

Perhaps you included the wrong samples or need to add an additional information to a PoolParty run that is currently underway. Since many modules run background processes, you will have to kill the entire script in this fashion: 

Enter the module that is running and determine processes it is associated with:  

> $ ps -aux | grep PPalign 

Kill all script processes:  

> $ killall PPalign

You may need to kill any additional lingering processes

> $ kill PID
