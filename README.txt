============================================
Dodecahedral L-spaces and hyperbolic 4-manifolds, Section 2
============================================

This is part of the code and data accompanying the paper

* Ludovico Battista, Leonardo Ferrari and Diego Santoro; Dodecahedral L-spaces and hyperbolic 4-manifolds, 2022.

It consists of several files written by Nathan Dunfield and downloaded by
https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LCYXPO 
(see [Dun] for more information) and some files written by Ludovico Battista.


checkLspace
=========
This is the only directory that contains original code. The other two directories
just contain the scripts that we run using code provided by Nathan Dunfield. Also
this directory contains some scripts and a database provided by Nathan Dunfield.

In this directory you can find the proofs that 6 manifolds contained in
CubicalOrientableClosedCensus(betti=0) are L-spaces, and the code that we used to prove that.
To have more information on how it works, see the paper and, in particular, its appendix.
To see the proofs, open the Jupyter worksheet 

checkLspace/CubicalOrientableClosedCensus/Notebook_proofs_rightangleddodman.ipynb


foliar
=========
This directory contains code written by Nathan Dunfield. We used this code to find
coorientable taut foliations on 23 manifolds contained in CubicalOrientableClosedCensus(betti=0),
thus proving that they are not L-spaces. We also found coorientable taut foliations with
vanishing euler class on 15 such manifolds, thus proving that they have orderable fundamental group
(see [Dun] for more informations). To see the scripts, open the Jupyter worksheet 
foliar/Taut_fol_van_eu_class.ipynb

non_orderability
=========
This directory contains code written by Nathan Dunfield. We used this code to prove that
6 manifolds contained in CubicalOrientableClosedCensus(betti=0), have non-orientable fundamental group.
To see the scripts, open the Jupyter worksheet 
non_orderability/check_proofs_non_orderability.ipynb


Running the programs
====================

The program requires many packages to work. The best way to install them all 
is by downloading the docker image computop/sage:8.6 (see [doc]) by giving the following commands:

>sudo docker pull computop/sage:8.6

We reccomend the version 8.6 of Sage because the code from [Dun] is written 
to run in this version (in particular, in Python2), and we also used this version.
Once installed the docker image, one can run it with: 

>sudo docker run -it -p 127.0.0.1:8888:8888 computop/sage:8.6

We reccomend to run the image in this way to have the possibility to use Jupyter after.
At this moment, one should get a bash like the following:

sage@[IDCONTAINER]:~$

We want to copy the code we have into the container. To do this, we open a new terminal,
we move in the directory where we downloaded the file "dod_manifolds.tar.xz" and we give
the following command:

>sudo docker cp ./dod_manifolds.tar.xz [IDCONTAINER]:home/sage/dod_manifolds.tar.xz

We now need to unzip the directory into the container. To do this, we go back to the first terminal and give

sage@[IDCONTAINER]:~$ tar -xf dod_manifolds.tar.xz

At this point we need to install the packages into the container.
To do this, give the following commands (in the terminal where the container is running):

cd checkLspace
sage -pip install .
cd ..
cd foliar
sage -pip install .
cd ..
cd non_orderability/quickdisorder
sage -pip install .
cd ..
cd ..


Now you installed 3 packages in the sage available in the container: quickdisorder,
foliar and checkLspace. To use them, launch sage by typing

sage

and give the following

import snappy, quickdisorder, foliar, checkLspace

In the code there are several .ipynb files that can
be opened through Jupyter and that contain  examples 
of how the code works. To open jupyter, give 

sage@[IDCONTAINER]:~$ sage -n jupyterlab

and then open http://localhost:8888/lab on your browser. You should be
able to move in the files in the container now. The .ipynb files that you can
find to see examples of the code and the proofs available in the paper are:

checkLspace/CubicalOrientableClosedCensus/Notebook_proofs_rightangleddodman.ipynb
foliar/Taut_fol_van_eu_class.ipynb
non_orderability/check_proofs_non_orderability.ipynb


If you want to go back to the container you used (and where you put the
downloaded files and installed packages), save the [IDCONTAINER], open a new terminal and run

sudo docker restart [IDCONTAINER]
sudo docker exec -it [IDCONTAINER] bash



Getting help
============

You are strongly encouraged to contact me:

  ludox73@gmail.com

if you are interested in anything you find here.

License
=======
The files written by Nathan Dunfield were released under the public domain, as per:
  https://creativecommons.org/publicdomain/zero/1.0/legalcode

As the author of the remaining parts, I (Ludovico Battista) 
hereby release this code and data into the public domain, as per:
  https://creativecommons.org/publicdomain/zero/1.0/legalcode
  
 Bibliography
 ==============
 
 [Dun] Nathan M. Dunfield, Floer homology, group orderability, and taut
  foliations of hyperbolic 3-manifolds, 2019
  
 [doc]  https://snappy.math.uic.edu/installing.html#kitchen-sink
        https://hub.docker.com/r/computop/sage/
