# LandmarkPicker
This Python script is used to manually select matching landmark pairs between two H&E stained images from 10x Visium tissue sections. They are saved as CSV files (y, x) and are designed to be further used for affine transformation in STalignCode.py. 

## Description
This code is a command line Python script that opens two H&E images, one source tissue image and one reference tissue image, and allows the user to click matching landmark points between them. For each landmark pair, the script first displays the source image and asks the user to click a landmark. It then displays the reference image and asks the user to click the corresponding matching landmark. This is repeated for the prior selected number of landmark pairs (recommended 6 to 10).

The output files are saved in a directory named after the two samples and can be passed directly to STalignCode.py using the --points1 and --points2 arguments.

landmark reccs: landmark pairs should be matching clear anatomical or structural points of the tissues visible in both H&E images. Such might include tissue corners, adventitia on the outside of the section, folds, holes and indentations, internal structures (lumen). They should ideally be spread out throughout the whole tissue section and not clustered in one area. 

## Dependencies 
Python 3.9 and higher
The following Python packages: pandas, matplotlib (pyplot and image), argparse, pathlib
LandmarkPicker.py is a self-contained script and can be run by downloading it and passing arguments through the command line, replacing the example paths with file locations:

python LandmarkPicker.py \
  --image1 /path/to/source_sample/spatial/tissue_hires_image.png \
  --image2 /path/to/reference_sample/spatial/tissue_hires_image.png \
  --sample_aligned SourceSample \
  --sample_reference ReferenceSample

After running the scripts asks the user to select the number of landmark pairs (6 to 10 is recommended), user should input the desired number of pairs. Afterwards the source image opens and the user is left to select the first of the pair, with the reference opening afterwards allowing the user to select the corresponding landmark. Repeat until all landmarks are selected. 

the code thus requires four arguments, those being: image1 and 2 (.png source and reference hires images), sample_aligned (name of the source sample being aligned), sample_reference (name of the reference sample)

optional arguments include: project_dir (main project directory where the landmark output folder will be created)

## Outputs 
Outputs are saved in a directory named after the two samples. The directory contains source landmark CSV, reference landmark CSV and combined landmark CSV for QC. For STalignCode.py only the first two are essential, the last one serves as a check whether landmark order has been preserved.

## Additional Notes and QC
-- matplotlib records clicked coordinates as x,y, but this script saves them as y,x because that is the coordinate format expected by STalignCode.py
-- for H&E landmark picking from tissue_hires_image.png, the resulting landmark coordinates should be used with: --coord_scale hires in STalignCode.py
-- in STalignCode.py would be passed as arguments pos1 (source tissue CSV) and pos2 (reference tissue CSV)