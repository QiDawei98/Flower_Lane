### WARNING: A bug has been reported by Liu Guoliang that occurs when the minimum value of sl2rho is greater than -0.4 a.u. A fix is in progress. ###
# Flower_Lane

A python script to perform Non Covalent Interactions index (NCI) analysis using electron density calculated by VASP (from CHGCAR files).  

Please cite:  [https://github.com/QiDawei98/Flower_Lane](https://github.com/QiDawei98/Flower_Lane)

## How to use?
Place your CHGCAR in the path specified in script then run it. The calculation may take 15 minutes or longer depending on your system.  

Open the resulting RDG (Reduced Density Gradient) file in Vesta, then use edit > edit data > volumetric data to load sl2rho file to surface coloring. Set the isosurface value to 0.5. Vesta will then render regions for non-covalent interactions, using color to indicate their nature (attractive, repulsive).

This project currently lacks documentation. I will add documentation and fix bugs as my schedule permits. In the meantime, if you encounter any difficulties, please feel free to email me with your CHGCAR file. 

## Algorithm
The program uses finite difference method on RBF interpolated data points. Visually, the results are indistinguishable from direct finite difference gradients calculated on cubic crystal systems using numpy.gradient and without interpolation.

Existing NCI implementations produce graphs that are visually less appealing than ours and do not accurately match results calculated directly from cubic systems. Furthermore, existing methods are known to cause unintuitive repulsion regions between ionic bonds.

Despite its simplicity, we believe our algorithm is superior to other currently available options.

## Acknowledgements
This is a collaborative work with Wang Xiu. We created the first fully functional prototype program together.

## Contact
Email me for **ANY** inquiry:  
[qidawei98@outlook.com](mailto:qidawei98@outlook.com)


