# FAQ

 __*Q:* The following error is printed upon program execution: ''Unable to generate or load fasta sequence file. Please check if 'bedtools' software is installed.'__  
 *A:* It seems the problem is related to missing bedtools software in the OS environment. In case of Ubuntu Linux, you need to execute the following command: `sudo apt-get install bedtools` which will install bedtools.


*⠀*



 __*Q:* The following error is printed upon program execution: 'Unable to convert the pdf file to png.'. No png (raster) image was generated.__
 
 *A:* You need to edit ImageMagick image processing permissions in your OS/user account. In Ubuntu Linux they should be at ‘/etc/ImageMagick-6/policy.xml’. You need to edit the file with root privileges (e.g. sudo) (using for example `sudo nano /etc/ImageMagick-6/policy.xml`). The line `<policy domain = "coder" rights = "none" pattern="PDF" />` should be changed by replacing `'rights = “none”'` with `‘rights = “read”'`.


*⠀*

 __*Q:* Something weird is happening, e.g. program is crashing with my data. What should I do?.__

 *A:* Please add `-debug` to the command line and post the details in [Issues](https://bitbucket.org/artegorov/svist4get/issues) or e-mail directly to _artyom**dot**egorov**AT**belozersky**dot**msu**dot**ru_.


### [Main page](https://github.com/art-egorov/svist4get)


#### [Quickstart guide](./QSGUIDE.md)

#### [Command-line parameters](/PARAMETERS.md)

#### [Configuration file parameters](./CONFIG.md)

#### [API usage examples](./API.md)

#### [Version log](./VERSION.md)
