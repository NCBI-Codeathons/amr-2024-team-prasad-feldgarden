# To install vadr follow this: https://github.com/ncbi/vadr/blob/master/documentation/install.md

## Installation notes for clean ubuntu 22.04 VM

```
curl -O https://raw.githubusercontent.com/ncbi/vadr/master/vadr-install.sh

# Instal unzip
sudo apt-get install unzip
# Instal autoconf
sudo apt-get install autoconf
# Instal cpanm
sudo apt-get install cpanminus
# instal Inline.pm
sudo cpanm install Inline
sudo cpanm install Inline::C
sudo apt install zlib1g
sudo cpanm install LWP::Simple
sudo cpanm install LWP::Protocol::https Mozilla::CA

bash vadr-install.sh linux
```

* Build you model using v-build
  * Choose an accession that will be your reference model (it can be a refseq or another reputable sequence)
  * For example https://www.ncbi.nlm.nih.gov/nuccore/HM000041, Klebsiella pneumoniae strain KP398/08
# Using v-build, build your model, it should take only a couple of minutes 
(information about v-build is here: https://github.com/ncbi/vadr/blob/master/documentation/build.md)

       $v-build.pl HM000041 HM000041

 # Then use v-annotate to compare related sequences to the model that you build
 (information about v-annotate  commands are here: https://github.com/ncbi/vadr/blob/master/documentation/advbuild.md)
    # First combine all files from the model that you built to one file (kp-models1)
    
    # If you have only one model sequence, this step isn't necessary, but I like to anyway
    
      # concatenate .minfo, .cm .fa and .hmm files:                                                                                                                                                                                    
        $ cat HM000041/*.vadr.minfo > kp-models1/kp.minfo
        $ cat HM000041/*.vadr.cm > kp-models1/kp.cm
        $ cat HM000041/*.vadr.fa > kp-models1/kp.fa
        $ cat HM000041/*.vadr.protein.hmm > kp-models1/kp.hmm
        
      # copy the blastdb files:
        $ cp HM000041/*.vadr.protein.fa* kp-models1/
        
      # prepare the library files:
        $ $VADREASELDIR/esl-sfetch --index kp-models1/kp.fa
        $ $VADRINFERNALDIR/cmpress kp-models1/kp.cm
        $ $VADRHMMERDIR/hmmpress kp-models1/kp.hmm
        $ $VADRBLASTDIR/makeblastdb -dbtype nucl -in kp-models1/kp.fa

  Run v-annotate.pl, and place those files into a new folder titled va-kp1 or (whatever you want)
  
        $ v-annotate.pl -f --out_stk --mdir kp-models1 --mkey kp kp-models1/kp.fa va-kp1

  The files that will give you the most information are the vadr.alt and the align.stk files in particular. 
