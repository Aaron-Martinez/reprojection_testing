# reprojection_testing
Old testing sandbox code for private reprojection project. Warning it is very messy

I created reprojection while I was working for the UGA Small Satellite Research Lab in college. The purpose of reprojection was to test an idea for a sub-algoritm to potentially be used in the structure from motion pipeline in SSRL's MOCI project. The MOCI small satellite needed to be able to perform structure from motion during low-earth orbit on an onboard GPU. Structure from motion is a technique to reconstruct a 3D scene from a series of two dimensional pictures taken from a camera at various locations or along a defined movement path.

This code contains some two-view examples of point clouds which are projected onto two cameras. They then should be "reprojected" back to the original 3D coordinates. This problem in the two-view case is referred to as triangulation. For N views a more advanced algorithm is needed such as performing least squares minimization of "reprojection error". 


This code makes use of the C++ numerical analysis library ALBLIB. See here: https://www.alglib.net/
