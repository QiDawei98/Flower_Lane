# Flower_Lane

A python script to perform Non Covalent Interactions index (NCI) analysis using electron density calculated by VASP (from CHGCAR files).

## How to use?
Place your CHGCAR in the path specified in script then run it. The calculation may take 15 minutes or longer depending on your system.  

Open the resulting RDG (Reduced Density Gradient) file in Vesta, then use edit > edit data > volumetric data to load sl2rho file to surce coloring. Set the isosurface value to 0.5. Vesta will then render regions for non-covalent interactions, using color to indicate their nature (attractive, repulsive).

Comments and documentation will be added soon, covering basic usage and details of the implementation.

## Algorithm
The program uses finite difference method on RBF interpolated data points. Visually, the results are indistinguishable from direct finite difference gradients calculated on cubic crystal systems using numpy.gradient and without interpolation.

Existing NCI implementations produce graphs that are visually less appealing than ours and do not accurately match results calculated directly from cubic systems. Furthermore, existing methods are known to cause unintuitive repulsion regions between ionic bonds.

Despite its simplicity, we believe our algorithm is superior to other currently available options.

## Acknowledgements
This is a collaborative work with Wang Xiu. We created the first fully functional prototype program together.

## Contact
Email me for **ANY** inquiry:  
[qidawei98@outlook.com](mailto:qidawei98@outlook.com)


