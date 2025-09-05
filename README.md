# MIL-GST
A Python library to compute mean intercept length (MIL) and gradient structure (GST) 2nd order tenors 

This projects computes 2nd order fabric tensors using either the mean intercept length method [1] or the gradient strucutre tensor [2]. These concepts are used, e.g., in tissue biomechanics [3-7] or  to quantify material microstructres for computational models and analyses.  

# Usage

## for calculatung MIL call
python3 microCT-to-microstructure.py -i image-data/ -b f -im example.mhd -t mil -w outputfolder/

OR

python3 microCT-to-microstructure.py -i image-data/ -b t -im example.mhd -t mil -w outputfolder/

## for calculating GST
python3 microCT-to-microstructure.py -i image-data/ -b f -im example.mhd -t gst -w outputfolder/

# References
[1] Harrigan, T. & Mann, R. Characterisation of Microstructural Anisotropy in Orthotropic Materials Using a Second Rank Tensor. Journal of Materials Research 19, 761–767 (1984).
[2] Tabor, Z. & Rokita, E. Quantifying anisotropy of trabecular bone from gray-level images. Bone 40, 966–972 (2007).
[3] Zysset, P. A Review of Morphology-Elasticity Relationships in Human Trabecular Bone: Theories and Experiments. Journal of Biomechanics 36, 1469–1485 (2003). 
[4] Wolfram, U., Schmitz, B., Heuer, F., Reinehr, M. & Wilke, H.-J. Vertebral trabecular main direction can be determined from clinical CT datasets using the gradient structure tensor and not the inertia tensor – a     case study. Journal of Biomechanics 42, 1390–1396 (2009).
[5] Wolfram, U., Gross, T., Pahr, D., Schwiedrzik, J., Wilke, H.-J. & Zysset, P. Fabric Based Tsai-Wu Yield Criteria for Vertebral Trabecular Bone in Stress and Strain Space. Journal of the Mechanical Behavior of       Biomedical Materials 15, 218–228 (2012).
[6] Schwiedrzik, J. & Zysset, P. An anisotropic elastic-viscoplastic damage model for bone tissue. Biomechanics and Modeling in Mechanobiology 12, 201–213 (2013).
[7] Schmidt J & Hartmaier A. A new texture descriptor for data-driven constitutive modeling of anisotropic plasticity. J Mater Sci 58, 14029–14050 (2023) doi: 10.1007/s10853-023-08852-2
