import numpy as np
import nibabel as nib
import pydicom as pyd
import argparse
from matplotlib import pyplot as plt
import matplotlib as mpl
import glob as glob
import os as os

# Interfile to NIfTI conversion tool for images reconstructed with CASToR
# See available options for conversion of 3D volumes, 3D+time, 3D+gating, 3D+time+gating
# Only some types of Interfile/NIfTI images are supported
# By default, use NIfTI version 1, as it is more widely supported
# For time frames images, as there is no actual dynamic support in NIfTI, write time information as embedded DICOM fields (used by Pmod), not very reliable

# input arguments
parser = argparse.ArgumentParser(description='Interfile to NIfTI conversion tool for images reconstructed with CASToR', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-fileBaseName', type=str, help='Name of the CASToR Interfile image without file extensions (.hdr/.img) and without frame/gate indications (e.g. _fr1), relative or full path', required=True)
parser.add_argument('-nbFrames', type=int, help='Number of time frames to convert: input = one Interfile file per frame, output = single NIfTI file containing all the frames, compatible with Pmod', default=1)
parser.add_argument('-nbRgates', type=int, help='Number of respiratory gates: input = one Interfile file per gate, output = one NIfTI file per gate ', default=1)
parser.add_argument('-nbCgates', type=int, help='Number of cardiac gates: input = one Interfile file per gate, output = one NIfTI file per gate', default=1)
parser.add_argument('-niftiVersion', type=int, help="NIfTI version, either 1 or 2, default is 1", default=1)
parser.add_argument('-convertFloatToInt', action='store_true', help="Convert image voxel values from floating point precision to uint16, which implies using rescale slope/intercept. The range of values in float should not be very different from the final range in uint16.", default=False)
parser.add_argument('-outFolder', type=str, help="Output folder, by default equal to the location of input files", default='')


args = parser.parse_args()
nbFrames = args.nbFrames
fileBaseName = args.fileBaseName
nbRgates = args.nbRgates
nbCgates = args.nbCgates
niftiVersion = args.niftiVersion
convertFloatToInt = args.convertFloatToInt
outFolder = args.outFolder

if niftiVersion!=1 and niftiVersion!=2:
  print("NIfTI version either 1 or 2")
  quit()

# loop on respiratory and cardiac gates if any
for rg in range(1,nbRgates+1):
  for cg in range(1,nbCgates+1):

    # automatic CASToR file name suffix in presence of gating
    suffix = ''
    if nbRgates>1:
      suffix += '_rg{:d}'.format(rg)
    if nbCgates>1:
      suffix += '_cg{:d}'.format(cg)

    # if dynamic image series (time frames), initialize additional variables and 
    # read all header information except the framing information from the first frame
    if nbFrames>1:
      # information on frame start times and frame duration for all frames
      frameStartTime = []
      frameDuration = []

      # paths to Interfile header and binary image for the first frame
      hdrPath = fileBaseName+'_fr1'+suffix+'.hdr'
      imgPath = fileBaseName+'_fr1'+suffix+'.img'
    else:
      # paths to Interfile header and binary image
      hdrPath = fileBaseName+suffix+'.hdr'
      imgPath = fileBaseName+suffix+'.img'

    # parse the Interfile header into a dictionary
    header = dict()
    with open (hdrPath) as f:
      l = f.readline()
      while l:
        l = f.readline()
        # split keys from values
        s = l.split(':=')
        if len(s)>1:
          # store key-value pair, trim any white spaces
          header[s[0].strip().replace('!','')]=s[1].strip()

    # check some fields, only simple 3D volume Interfile files are supported currently
    if 'number of frame groups' in header.keys() and int(header['number of frame groups'])>1:
      print ('Number of frame groups > 1, not supported yet!')
      quit()

    # number format for voxel values
    if header['number format']=='short float' and header['number of bytes per pixel']=='4':
      precision = np.float32
    elif header['number format']=='long float' and header['number of bytes per pixel']=='8':
      precision = np.double
    else:
      print ('Number format not supported, only floating point number formats are supported!')

    # image dimensions (X Y Z)
    nbDim = int(header['number of dimensions'])
    # number of voxels for each dimension
    dimNbVox = np.zeros(nbDim, np.int32)
    # voxel size in mm for each dimension
    dimVoxSize = np.zeros(nbDim, np.float32)
    # offset for the field-of-view in mm
    offset = np.zeros(nbDim, np.float32)
    for d in range(0,nbDim):
      dimNbVox[d]=int(header['matrix size [{:d}]'.format(d+1)])
      dimVoxSize[d]=float(header['scaling factor (mm/pixel) [{:d}]'.format(d+1)])
      offset[d]=float(header['first pixel offset (mm) [{:d}]'.format(d+1)])

    # bed offset in Z, if available (if available in the CASToR datafile used for reconstruction, it should correspond to the DICOM positioning fields in the raw data)
    # TODO not sure this always works correctly
    if 'horizontal bed relative position (mm)' in header.keys():
      offset[2] -= float(header['horizontal bed relative position (mm)'])

    # Affine transform for NIfTI
    affine=np.array([[-dimVoxSize[0],0,0,0.5*(dimNbVox[0]-1)*dimVoxSize[0]+offset[0]],
    [0,dimVoxSize[1],0,-0.5*(dimNbVox[1]-1)*dimVoxSize[1]-offset[1]],
    [0,0,dimVoxSize[2],-0.5*(dimNbVox[2]-1)*dimVoxSize[2]-offset[2]],
    [0,0,0,1]])

    # read dynamic (time frame) information and images from other frames
    if nbFrames>1:
      # Actual image matrix (X Y Z T)
      im = np.zeros((dimNbVox[0], dimNbVox[1], dimNbVox[2], nbFrames),np.float32)  

      # load each frame and store the image matrix and the frame start time and duration 
      # the order of dimensions of the image matrix is inverted to match NIfTI
      for fr in range(0,nbFrames):
        # name of the frame files
        hdrPath = fileBaseName+'_fr{:d}'.format(fr+1)+suffix+'.hdr'
        imgPath = fileBaseName+'_fr{:d}'.format(fr+1)+suffix+'.img'

        # parse the Interfile header into a dictionary
        header = dict()
        with open (hdrPath) as f:
          l = f.readline()
          while l:
            l = f.readline()
            # split keys from values
            s = l.split(':=')
            if len(s)>1:
              # store key-value pair, trim any white spaces
              header[s[0].strip().replace('!','')]=s[1].strip()
        
        # load image matrix from the binary .img Interfile file and invert the order of dimensions to match NIfTI
        im[:,:,:,fr] = np.transpose(np.fromfile(imgPath, dtype=precision, count=-1).reshape(dimNbVox[::-1].tolist()))
        # frame start time information in s
        frameStartTime.append(float(header['image start time (sec)']))
        # frame start duration in ms, apparently required by Pmod
        frameDuration.append(float(header['image duration (sec)'])*1000.)
    else:
      # load image matrix from the binary .img Interfile file and invert the order of dimensions to match NIfTI
      im = np.transpose(np.fromfile(imgPath, dtype=precision, count=-1).reshape(dimNbVox[::-1].tolist()))


    # convert voxel values from float to uint16
    # compute rescale slope/intercept
    rescaleSlope=1.
    rescaleIntercept=0.
    if convertFloatToInt:
      print ("Initial range of image voxel values (format {}): min={:g} max={:g}".format(precision, np.min(im), np.max(im)))
      rescaleIntercept = np.min(im)
      im = im-rescaleIntercept

      rescaleSlope = np.max(im)/65535.
      if rescaleSlope>2.:
        print("Warning! The rescale slope is surprisingly high ({:g}), so the conversion from float to uint16 might compromise image voxel values. Check if your CASToR image has some very noisy values, at image edges for instance. If this is the case, apply a mask or similar to avoid very high unnecessary voxel values! ".format(rescaleSlope))
      im = im/rescaleSlope
      #convert image to uint16
      im = np.uint16(im) 
      print ("Final range of image voxel values (format uint16): min={:d} max={:d}".format(np.min(im), np.max(im)))
      
    # create the final NIfTI image
    if niftiVersion==1:
      niftiIm = nib.Nifti1Image(im, affine)
    elif niftiVersion==2:
      niftiIm = nib.Nifti2Image(im, affine)

    # set units
    niftiIm.header.set_xyzt_units('mm', 'sec')
    niftiIm.header.set_slope_inter(rescaleSlope,rescaleIntercept)
      
    # store information for the time dimension in DICOM extension, used by Pmod
    if nbFrames>1 and niftiVersion==1:    
      dicomds = pyd.Dataset()
      dicomds.add_new((0x0054,0x1001),'CS','Bq/ml')
      print("Warning! The unit for voxel values is set to Bq/ml by default, but we can't know if this is true from the current Interfile!")
      dicomds.add_new((0x0055,0x0010),'LO','PMOD_1')
      dicomds.add_new((0x0055,0x1001),'FD',frameStartTime)
      dicomds.add_new((0x0055,0x1004),'FD',frameDuration)
      dcmext = nib.nifti1.Nifti1DicomExtension(2, dicomds) # Use DICOM ecode 2
      niftiIm.header.extensions.append(dcmext)

    # save the final NIfTI image to .nii file
    if outFolder:
      inFolderName, inFileName = os.path.split(fileBaseName)
      nib.save(niftiIm, outFolder+inFileName+suffix+'.nii')
    else:
      nib.save(niftiIm, fileBaseName+suffix+'.nii')

