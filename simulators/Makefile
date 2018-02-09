# Downloads and compiles the required simulators.

SCRM_VERSION=1.6.1
MSMS_VERSION=3.2rc-b163
.PHONY: all clean
all: scrm ms msms.jar macs cosi2

cosi2:
	wget http://www.broadinstitute.org/mpg/cosi2/coalescent_linux_static
	mv coalescent_linux_static cosi2
	chmod u+x cosi2
	# NOTE for cosi2 we must set LC_ALL=C in the shell or it won't 
	# work.
	export LC_ALL=C

macs: macs-0.5d.tar.gz
	# Downoad for macs no longer works. Get from 
	# https://code.google.com/archive/p/macs/downloads
	# wget http://macs.googlecode.com/files/macs-0.5d.tar.gz
	tar -zxvf macs-0.5d.tar.gz
	mv macs macs_dir
	cd macs_dir && make
	cp macs_dir/macs ./ 

scrm:
	wget https://github.com/scrm/scrm/releases/download/v${SCRM_VERSION}/scrm-src.tar.gz
	tar -zxvf scrm-src.tar.gz
	cd scrm-${SCRM_VERSION} && ./configure
	make -C scrm-${SCRM_VERSION} 
	cp scrm-${SCRM_VERSION}/scrm ./

ms: ms.tar.gz
	# ms download doesn't work any more because the website won't give 
	# permanent links. Need to download from 
	# Http://home.uchicago.edu/rhudson1/source/mksamples.html
	# first then build.
	tar -zxvf ms.tar.gz
	cd msdir && gcc -O3 -o ms ms.c streec.c rand1.c -lm
	cp msdir/ms ./

msms.jar:
	wget http://www.mabs.at/ewing/msms/msms${MSMS_VERSION}.jar
	mv msms${MSMS_VERSION}.jar msms.jar

clean:
	rm -fR scrm*
	rm -fR ms*
	rm -fR macs*
	rm -fR cosi2*
