**Pynpre**
======

## Description

**Pynpre* is a python package for denoising and interpolation of multi-dimensional multi-channel seismic data. This package has a variety of applications in both exploration and earthquake seismology.

## Reference
Wang et al., 2021, Nonstationary predictive filtering for seismic random noise suppression - A tutorial, Geophysics, 86(3), W21–W30. 

Chen et al., 2021, 5D de-aliased seismic data interpolation using non-stationary prediction error filter, Geophysics, 86(5), V419–V429.

BibTeX:

	@article{npre,
	  title={Nonstationary predictive filtering for seismic random noise suppression - A tutorial},
	  author={Hang Wang and Wei Chen and Weilin Huang and Shaohuan Zu and Xingye Liu and Liuqing Yang and Yangkang Chen},
	  journal={Geophysics},
	  volume={86},
	  number={3},
	  issue={3},
	  pages={W21–W30},
	  doi={10.1190/geo2020-0368.1},
	  year={2021}
	}

	@article{npre5d,
	  title={5D de-aliased seismic data interpolation using non-stationary prediction error filter},
	  author={Yangkang Chen and Sergey Fomel and Hang Wang and Shaohuan Zu},
	  journal={Geophysics},
	  volume={86},
	  number={5},
	  issue={5},
	  pages={V419–V429},
	  doi={10.1190/geo2020.0540.1},
	  year={2021}
	}

-----------
## Copyright
    The pynpre developing team, 2021-present
-----------

## License
    MIT License   

-----------

## Install
Using the latest version (Suggested)

    git clone https://github.com/chenyk1990/pynpre
    cd pynpre
    pip install -v -e .
    
or using Pypi

    pip install pynpre

or
    pip install git+https://github.com/chenyk1990/pynpre
    
-----------
## Examples
    The "demo" directory contains all runable scripts to demonstrate different applications of pynpre. 

-----------
## Dependence Packages
* scipy 
* numpy 
* matplotlib

-----------
## Modules
    xxx.py  -> description
    
-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Yangkang Chen
    chenyk2016@gmail.com

-----------
## Gallery
The gallery figures of the pynpre package can be found at
    https://github.com/chenyk1990/gallery/tree/main/pynpre
Each figure in the gallery directory corresponds to a DEMO script in the "demo" directory with the exactly the same file name.

These gallery figures are also presented below. 

DEMO1 (test_pynpre_syn2d.py)

<img src='https://github.com/chenyk1990/gallery/blob/main/pynpre/test_pynpre_syn2d.png' alt='DEMO1' width=960/>

DEMO2 (test_pynpre_syn3d.py) This example is to show that there are no perfect denoising methods. One can work on a certain type of data, but not all. The plane-wave synthetic example is biased towards the damped rank-reduction (DRR) method. However, NPRE works perfectly on curving and non-stationary seismic data.

<img src='https://github.com/chenyk1990/gallery/blob/main/pynpre/test_pynpre_syn3d.png' alt='DEMO2' width=960/>

DEMO3 (test_pynpre_real3d.py) 

<img src='https://github.com/chenyk1990/gallery/blob/main/pynpre/test_pynpre_real3d.png' alt='DEMO3' width=960/>





# Below are other functions (Time-frequency analysis) enabled by the Pynpre package
Generated by [demos/test_pynpre_syn1d_stltft](https://github.com/chenyk1990/pynpre/tree/main/demos/test_pynpre_syn1d_stltft.py and https://github.com/chenyk1990/pyntfa/tree/main/demos/test_pyntfa_syn1d_stltft.py)

<img src='https://github.com/chenyk1990/gallery/blob/main/pynpre/test_pyntfa_syn1d_stltft_stft.png' alt='comp' width=960/>

<img src='https://github.com/chenyk1990/gallery/blob/main/pynpre/test_pyntfa_syn1d_stltft_st.png' alt='comp' width=960/>

<img src='https://github.com/chenyk1990/gallery/blob/main/pynpre/test_pyntfa_syn1d_stltft_ntfa.png' alt='comp' width=960/>


