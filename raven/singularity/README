************************
Using singularity images
************************

Setup
-----

* install singularity locally: https://www.sylabs.io/guides/2.6/user-guide/installation.html#install-on-linux
	* Perform an installation from source, option 1
	* version 2.5.2 is fine
	* Some external packages are required, such as libarchive-devel

Test using a Raven dataset
--------------------------

* In your working directory create a folder 'data' and fill it with the raven input data:
	./data/Irond.rvi
	./data/Irond.rvp
	...
	
* Create an output folder, e.g. 

	$ mkdir data_out
	
* Set an environment variable, i.e. (for bash):

	$ export SINGULARITY_NOHTTPS=true
	
* Invoke singularity as follows:

	$ singularity run --bind ./data:/data  --bind ./data_out:/data_out:rw shub://132.217.141.54/hydro/raven:latest <dataset_name>
	where <dataset_name> may be e.g. Irond.
	
A progress bar appears and indicates that the raven singularity image is being downloaded. When download is complete, the raven 
application executes and generates its results inside the data_out folder.


******************************
Maintaining singularity images
******************************

Get Raven source code (v2.8.1):
./get_raven.sh

Build singularity image and push it to registry (needs sudo):
./prepare_simg.sh

In order to be able to push images to the registry (located at http://132.217.141.54/), one must 
1- log in to the registry via github
2- go to the user menu (right to the menu 'Tools')
3- select the menu item 'Token' and follow the instructions
 
