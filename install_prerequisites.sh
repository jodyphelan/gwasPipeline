mkdir bin
# Plink
wget --no-check-certificate https://www.cog-genomics.org/static/bin/plink170906/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
mv plink bin/
rm toy* prettify

# Beagle
wget --no-check-certificate https://faculty.washington.edu/browning/beagle/beagle.08Jun17.d8b.jar
mv beagle.08Jun17.d8b.jar bin/

# htslib
cd htslib
autoheader
autoconf
./configure --disable-lzma;make
mv tabix bgzip ../bin/
cd ../
