#!/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 14:19:43 2025

@author: uw62
"""

'''
@ Author: Uwe Wolfram

forked of code from Marta Peña Fernández

------------------------------------------------------------------------------
The script reads a .mhd image, either in TIFF or binary format, and generates
a input file for simple compression or tension testing in Abaqus
'''


import SimpleITK as sitk
import os
from os import path
from glob import glob
from datetime import date

import numpy as np

#from distutils.util import strtobool
from str2bool import str2bool
import argparse

from skimage.measure import label   
from skimage.measure import marching_cubes
from scipy.ndimage import gaussian_filter
import pyvista as pv   # pip install pyvista

def getInputs():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--inp',
                        default='./', # i.e. local folder
                        type=str,
                        help='input image folder'
                        )   
    parser.add_argument('-im', '--image',
                        default='Trabecular-Bone-Test.mhd',
                        type=str,
                        help='image used for simulation'
                        )
    parser.add_argument('-b', '--binary',
                        default=str(False),
                        type=str2bool,
                        help='binary image, i.e. consisting of [o,1]'
                        )
    parser.add_argument('-t', '--tensor',
                        default='mil',
                        type=str,
                        help='choice is either mil or gst'
                        )
    parser.add_argument('-w', '--write',
                        default='./', # i.e. local folder
                        type=str,
                        help='folder to write abaqus input file'
                        )

    args = parser.parse_args()
    inp = args.inp
    binary = args.binary
    write = args.write
    im = args.image
    tensor = args.tensor

    return inp, write, binary, im, tensor

def getLargestCC(segmentation):

    labels = label(segmentation, background=False)
    unique, counts = np.unique(labels, return_counts=True)
    list_seg=list(zip(unique, counts))[1:] # the 0 label is by default background so take the rest
    largest=max(list_seg, key=lambda x:x[1])[0]
    labels_max=(labels == largest).astype(int)
    
    return labels_max

def mhdReader(mhdFile):
    
    itkimage = sitk.ReadImage(mhdFile)    
    V = np.array(sitk.GetArrayFromImage(itkimage),dtype=bool)    
    spacing = np.array(list(reversed(itkimage.GetSpacing())))    
    V = np.transpose(V, (2, 1, 0))
    
    return V, spacing

def threshold(mhdFile):
    
    # read MHD
    itkimage = sitk.ReadImage(mhdFile)  
    
    # apply median filter, radius=2
    medianFilter = sitk.MedianImageFilter()
    medianFilter.SetRadius(2)
    medianimage  = medianFilter.Execute(itkimage)
    
    # apply otsu filter
    otsu_filter = sitk.OtsuThresholdImageFilter()
    otsu_filter.SetInsideValue(0)
    otsu_filter.SetOutsideValue(1)
    
    image = otsu_filter.Execute(medianimage)
    
    for i in range(0,2):
        erode_filter = sitk.BinaryErodeImageFilter()
        image = erode_filter.Execute(image)
        
        dilate_filter = sitk.BinaryDilateImageFilter()
        image = dilate_filter.Execute(image)
        
        dilate_filter = sitk.BinaryDilateImageFilter()
        image = dilate_filter.Execute(image)
        
        erode_filter = sitk.BinaryErodeImageFilter()
        image = erode_filter.Execute(image)
    
    # get 3D array of image - if connectivity filter is NOT used
    V = sitk.GetArrayFromImage(image)
    
    V = np.transpose(V, (2, 1, 0))
    
    # read spacing from MHD
    spacing = np.array(list(reversed(itkimage.GetSpacing())))
    
    # extract original gray scale image as an array as well
    I = sitk.GetArrayFromImage(itkimage)
    
    I = np.transpose(I, (2, 1, 0))
    
    return V, I, spacing

def compute_bvtv(V):
    
    tv = V.shape[0]*V.shape[1]*V.shape[2]
    bv = (V == True).sum()
            
    return bv/tv



def _spherical_fibonacci_points(n):
    """
    Nearly uniform points on S^2 (unit sphere) using spherical Fibonacci lattice.
    Returns n x 3 array of unit vectors.
    """
    # Golden ratio
    phi = (1 + 5 ** 0.5) / 2
    i = np.arange(n) + 0.5
    theta = 2 * np.pi * i / phi
    z = 1 - 2 * i / n
    r = np.sqrt(np.maximum(0.0, 1 - z * z))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    dirs = np.stack([x, y, z], axis=1)
    # Remove antipodal duplicates statistically by enforcing u_z >= 0 then mirroring in fit
    return dirs / np.linalg.norm(dirs, axis=1, keepdims=True)

def compute_mil_tensor(
    volume: np.ndarray,
    spacing=(1.0, 1.0, 1.0),
    n_directions: int = 600,
    iso_level: float = 0.5,
    ):
    """
    Compute the 3x3 Mean Intercept Length (MIL) tensor from a segmented 3D volume.

    Parameters
    ----------
    volume : np.ndarray
        3D numpy array (bool or 0/1) with True/1 = solid (foreground), False/0 = void (background).
    spacing : tuple of float
        Physical voxel spacing (dz, dy, dx). Use (1,1,1) if unknown.
    n_directions : int
        Number of approximately-uniform sphere directions to sample (>= 200 recommended).
    iso_level : float
        Isosurface level for marching cubes (0.5 is standard for binary arrays).

    Returns
    -------
    M : (3,3) np.ndarray
        The MIL tensor (symmetric positive-definite). Eigenvectors give principal fabric directions;
        eigenvalues relate to principal MILs.
        
        evals is length-3 array of eigenvalues.
        evecs is a matrix with eigenvectors as columns, not rows.
    info : dict
        Extra diagnostics (principal MILs, eigenvectors, raw directional samples).
    """
    # --- Validations ---
    if volume.ndim != 3:
        raise ValueError("volume must be a 3D array")
    vol_bin = (volume > 0).astype(np.uint8)



    # Physical voxel volume and total physical volume
    dz, dy, dx = spacing
    voxel_vol = float(dz) * float(dy) * float(dx)
    V_total = vol_bin.size * voxel_vol
    V_solid = vol_bin.sum() * voxel_vol
    if V_solid == 0 or V_solid == V_total:
        raise ValueError("Volume must contain both solid and void to define an interface.")


    # --- Extract interface via marching cubes (solid-void surface) ---
    # This produces a triangle mesh approximating the isosurface at level iso_level
    # try:
    #     from skimage.measure import marching_cubes
    # except Exception as e:
    #     raise ImportError(
    #         "This function requires scikit-image (skimage). Please install it (e.g., pip install scikit-image)."
    #     ) from e

   
    # marching_cubes expects array in (z, y, x) with spacing=(dz, dy, dx)
    verts, faces, _, _ = marching_cubes(vol_bin, level=iso_level, spacing=spacing)

    # --- Compute per-face unit normals and areas in physical units ---
    # Face vertices
    v0 = verts[faces[:, 0]]
    v1 = verts[faces[:, 1]]
    v2 = verts[faces[:, 2]]
    # Edge vectors
    e1 = v1 - v0
    e2 = v2 - v0
    # Unnormalized normals = cross product; area = 0.5 * ||cross||
    cross = np.cross(e1, e2)
    face_area = 0.5 * np.linalg.norm(cross, axis=1)  # shape (F,)
    # Unit normals; avoid divide-by-zero for degenerate triangles
    with np.errstate(invalid='ignore', divide='ignore'):
        n_face = cross / (2.0 * face_area[:, None])
    # Replace degenerate normals with zeros
    n_face[~np.isfinite(n_face)] = 0.0

    # --- Sample directions on the unit sphere ---
    U = _spherical_fibonacci_points(n_directions)  # (K, 3)

    # --- Directional surface-area density S_dir(u) per unit volume ---
    # S_dir(u) ≈ (1/V_total) * sum_faces a_f * |n_f · u|
    # Vectorized: for each u, compute |N u| via dot
    # N: (F,3), U: (K,3) -> dots: (F,K)
    dots = n_face @ U.T                      # (F, K)
    S_dir = (face_area[:, None] * np.abs(dots)).sum(axis=0) / V_total  # (K,)

    # Avoid zeros (perfectly flat cases numerically)
    eps = np.finfo(float).eps
    S_dir = np.maximum(S_dir, 10 * eps)

    # --- Directional MIL values ---
    MIL_dir = (2.0 * V_solid) / S_dir    # (K,)

    # --- Fit quadratic form: 1 / MIL(u) ≈ u^T B u (B symmetric 3x3) ---
    # Build design matrix for each u = (ux,uy,uz): [ux^2, uy^2, uz^2, 2uxuy, 2uxuz, 2uyuz]
    ux, uy, uz = U[:, 0], U[:, 1], U[:, 2]
    X = np.column_stack([
        ux * ux,
        uy * uy,
        uz * uz,
        2.0 * ux * uy,
        2.0 * ux * uz,
        2.0 * uy * uz,
    ])  # shape (K, 6)

    y = 1.0 / MIL_dir  # target

    # Least squares solve
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)  # (6,)

    # Reconstruct symmetric B
    B = np.array([
        [beta[0], beta[3] / 2.0, beta[4] / 2.0],
        [beta[3] / 2.0, beta[1], beta[5] / 2.0],
        [beta[4] / 2.0, beta[5] / 2.0, beta[2]],
    ])

    # Regularize to ensure positive-definite (tiny jitter if needed)
    # (Should already be PD in well-behaved data)
    # same as my old one
    w, V = np.linalg.eigh(B)
    w = np.maximum(w, 1e-12)
    B_pd = (V * w) @ V.T

    # MIL tensor is the inverse of the fitted quadratic form
    M = np.linalg.inv(B_pd)

    # Principal MILs and directions (eigendecomposition of M)
    evals, evecs = np.linalg.eigh(M)
    # Sort ascending by eigenvalue (largest MIL associated with material main axis)
    # normalisation sum(evals) = 3
    idx = np.argsort(evals)[::1]
    evals = evals[idx]
    evecs = evecs[:, idx]

    # norming evals so that tr(G) = 3    
    evals = 3*evals/np.sum(evals)
    
    
    # reconstruct fabric tensor. Note: evecs[:, 0] is the first eigenvector. evecs[0] is just the 
    # first row of the whole eigenvector matrix. This is due to the shape convention that np.linalg.eigh uses 
    # for eigenvectors
    M = evals[0]*np.outer(evecs[:,0], evecs[:,0]) + evals[1]*np.outer(evecs[:,1], evecs[:,1]) + evals[2]*np.outer(evecs[:,2], evecs[:,2])
    
    info = dict(
        principal_MILs=evals,
        principal_dirs=evecs,           # columns correspond to principal_MILs
        sampled_dirs=U,
        MIL_per_dir=MIL_dir,
        B=B_pd,
        solid_fraction=float(V_solid / V_total),
        total_surface_area=float(face_area.sum()),
    )
    
    
    return M, info
 

def compute_gst_tensor(
    image: np.ndarray,
    spacing=(1.0, 1.0, 1.0),
    sigma_e=1.0,
    sigma_i=2.0
    ):
    """
    Compute the 3x3 gradient structure tensor (GST) from a gray valued 3D volume.

    Parameters
    ----------
    volume : np.ndarray
        3D numpy array with gray scales.
    spacing : tuple of float
        Physical voxel spacing (dz, dy, dx). Use (1,1,1) if unknown.
    sigma_e : float
        Derivative scale (Gaussian smoothing before gradients), in *physical units*.
        How much you smooth before taking gradients (noise control).
    sigma_i : float
        Integration scale (smoothing after outer-product), in *physical units*.
        How much you smooth the per-voxel tensors after outer-product (spatial coherence).
        
    Notes
    -----
    Choosing σ’s:
        Start with sigma_e ≈ 1–2 voxels and sigma_i ≈ 2–4 voxels (in physical units).
        Increase sigma_e if your data is noisy; increase sigma_i if you want orientations to vary smoothly across space.
    Physical correctness: 
        Gradients are computed with np.gradient(..., dz, dy, dx) so derivatives are in per-unit-length; σ’s are specified 
        in physical units and converted internally to pixel units per axis, so the result is spacing-aware.
    Binary inputs: 
        For segmented images, a mild sigma_e helps avoid staircase gradients; the tensor then forms on the interfaces, 
        which is exactly what we want.
    
    Returns
    -------
    G : (3,3) np.ndarray
        The MIL tensor (symmetric positive-definite). Eigenvectors give principal fabric directions;
        eigenvalues relate to principal MILs.
        
        evals is length-3 array of eigenvalues.
        evecs is a matrix with eigenvectors as columns, not rows.
        
    info : dict
        Extra diagnostics (principal MILs, eigenvectors, raw directional samples).
    """
    # --- Validations ---
    if image.ndim != 3:
        raise ValueError("image must be a 3D array")
    vol_bin = (image > 0).astype(np.uint8)
    
    I = np.asarray(image, dtype=np.float32)
    dz, dy, dx = spacing
    
    # Convert σ from physical units to pixel units per axis
    sigma_e_pix = (sigma_e / dz, sigma_e / dy, sigma_e / dx)
    sigma_i_pix = (sigma_i / dz, sigma_i / dy, sigma_i / dx)

    
    # smooth at derivative scale
    I_s = gaussian_filter(I, sigma=sigma_e_pix, mode="nearest")
    
    # gradients in *physical units*
    Gz, Gy, Gx = np.gradient(I_s, dz, dy, dx, edge_order=2)
    
    # calculating tensor components
    gzz=Gz*Gz
    gyy=Gy*Gy
    gxx=Gx*Gx
    
    gyz=Gy*Gz
    gzx=Gz*Gy
    gxy=Gx*Gy
    
    # integrating/smooth at integration scale ---
    gzz = gaussian_filter(gzz, sigma=sigma_i_pix, mode="nearest")
    gyy = gaussian_filter(gyy, sigma=sigma_i_pix, mode="nearest")
    gxx = gaussian_filter(gxx, sigma=sigma_i_pix, mode="nearest")
    gyz = gaussian_filter(gyz, sigma=sigma_i_pix, mode="nearest")
    gzx = gaussian_filter(gzx, sigma=sigma_i_pix, mode="nearest")
    gxy = gaussian_filter(gxy, sigma=sigma_i_pix, mode="nearest")
    
    G = np.array([[gxx.mean(), gxy.mean(), gzx.mean()],
                 [gxy.mean(), gyy.mean(), gyz.mean()],
                 [gzx.mean(), gyz.mean(), gzz.mean()]])
    
    # Principal MILs and directions (eigendecomposition of M)
    evals, evecs = np.linalg.eigh(G)

    # norming evals so that tr(G) = 3    
    evals = 3*evals/np.sum(evals)
    
    
    # reconstruct fabric tensor. Note: evecs[:, 0] is the first eigenvector. evecs[0] is just the 
    # first row of the whole eigenvector matrix. This is due to the shape convention that np.linalg.eigh uses 
    # for eigenvectors
    G = evals[0]*np.outer(evecs[:,0], evecs[:,0]) + evals[1]*np.outer(evecs[:,1], evecs[:,1]) + evals[2]*np.outer(evecs[:,2], evecs[:,2])
        
    """
    # Sort ascending by eigenvalue (largest MIL associated with material main axis)
    # normalisation sum(evals) = 3
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    """
    
    
    """
    
    CHECK REMOVAL OF THE COMMENTED ABOVE
    
    """
    
    evals = 3*evals/np.sum(evals)
    
    info = dict(
        principal_GSTs=evals,
        principal_dirs=evecs,           # columns correspond to principal_MILs
        )
    
    return G, info
    

def tensor_to_vtk(tensor: np.ndarray,
                  filename="tensor_ellipsoid.vtk",
                  n_theta=60,
                  n_phi=30,
                  scale=1.0):
    """
    Convert a 3x3 symmetric tensor into a 3D ellipsoid mesh and save as VTK.

    Parameters
    ----------
    tensor : (3,3) ndarray
        Symmetric 2nd order tensor (e.g. MIL or GST).
    filename : str
        Output VTK filename.
    n_theta : int
        Number of divisions around azimuth.
    n_phi : int
        Number of divisions in polar angle.
    """

    # 1) Eigen-decomposition
    evals, evecs = np.linalg.eigh(tensor)
    idx = np.argsort(evals)[::-1]
    evals, evecs = evals[idx], evecs[:, idx]
    
    # Protect against negatives (numerical noise)
    evals = np.clip(evals, 1e-12, None)

    # 2) Parametric sphere
    theta = np.linspace(0, 2*np.pi, n_theta)
    phi = np.linspace(0, np.pi, n_phi)
    theta, phi = np.meshgrid(theta, phi)

    xs = np.sin(phi) * np.cos(theta)
    ys = np.sin(phi) * np.sin(theta)
    zs = np.cos(phi)
    sphere = np.stack([xs, ys, zs], axis=-1)
    
    # 3) Scale along eigen-directions
    axes_lengths = scale * np.sqrt(evals)   # semi-axes
    ellipsoid = sphere @ np.diag(axes_lengths) @ evecs.T
    
    X, Y, Z = ellipsoid[..., 0], ellipsoid[..., 1], ellipsoid[..., 2]
    
    # 3) Build ellipsoid surface ---
    surf = pv.StructuredGrid(X, Y, Z)

    # 4) Build principal axes arrows 
    origin = np.zeros(3)
    arrows = []
    colors = ["red", "green", "blue"]
    for i in range(3):
        vec = evecs[:, i] * axes_lengths[i]
        arrow = pv.Arrow(start=origin, direction=vec, scale="auto")
        arrow["color"] = np.full(arrow.n_points, i)  # encode index as scalar
        arrows.append(arrow)

    # 5) Combine into MultiBlock dataset ---
    multiblock = pv.MultiBlock()
    multiblock["Ellipsoid"] = surf
    for i, arrow in enumerate(arrows):
        multiblock[f"Axis_{i}"] = arrow

    # 6) Save everything into one file ---
    multiblock.save(filename)
    print(f"Ellipsoid + axes written to {filename}")

    return multiblock


    """
    # alternatiove version where M is written as ODF r = sqrt(x^T M x) and not ellipsoid
    # 1) Eigen-decomposition
    evals, evecs = np.linalg.eigh(tensor)
    idx = np.argsort(evals)[::-1]
    evals, evecs = evals[idx], evecs[:, idx]

    # Protect against negatives (numerical noise)
    evals = np.clip(evals, 1e-12, None)

    # 2) Parametric sphere
    theta = np.linspace(0, 2*np.pi, n_theta)
    phi = np.linspace(0, np.pi, n_phi)
    theta, phi = np.meshgrid(theta, phi)

    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    dirs = np.stack([x, y, z], axis=-1)

    # 3) Quadratic form scaling
    q = np.einsum('...i,ij,...j->...', dirs, tensor, dirs)
    r = np.sqrt(np.clip(q, 0, None))

    X = r * x
    Y = r * y
    Z = r * z

    # 4) Rotate into eigenbasis (align ellipsoid)
    R = evecs
    XYZ = np.stack([X, Y, Z], axis=-1)
    XYZ_rot = XYZ @ R.T
    Xr, Yr, Zr = XYZ_rot[...,0], XYZ_rot[...,1], XYZ_rot[...,2]

    # 5) Build pyvista mesh
    surf = pv.StructuredGrid(Xr, Yr, Zr)

    # 6) Save to VTK
    surf.save(filename)
    print(f"Ellipsoid written to {filename}")

    return surf
    """

def unique(Nodes):
    
    # function for returning repeated indices of unique nodes
    
    #label all duplicate nodes with ids
    unique,index,inv = np.unique(Nodes[:,1:4], return_index=True,
                                 return_inverse=True, axis=0)
    
    UniqueNodesMapped = np.zeros((np.size(index,0),4))
    UniqueNodesMapped[:,0] = np.arange(1,np.size(index)+1)
    UniqueNodesMapped[:,1:4] = Nodes[np.sort(index),1:4]
    
    ids = np.arange(len(unique)) #create node ids
    ids = ids[:,np.newaxis] #make it 2d array for concatenation with 'unique'
    
    unique_ids = np.concatenate((ids, unique), axis=1)
    NodesMapped = unique_ids[inv]
    
    #find first instances of each unique node
    unique,index = np.unique(Nodes[:,1:4], return_index=True, axis=0)
    
    first_instance = np.sort(index)
    
    #create an array of unique node ids and corresponding indices
    node_ids = NodesMapped[first_instance,0]
    node_ids = node_ids[:,np.newaxis]
    node_ids = np.concatenate((ids,node_ids), axis=1)
    node_ids = node_ids[node_ids[:,1].argsort()]
    
    #use indices to map unique node ids to nodes
    
    new_ids = node_ids[NodesMapped[:,0]]
    
    NodesMapped[:,0] = new_ids[:,0]

    return NodesMapped, UniqueNodesMapped   

def cubicMesh(Vol,spacing):
    
    
    V = getLargestCC(Vol)
    
       
    # Connected components only 
    
    
    # reshape V as nx3 array
    PtClouds = np.array(np.where(V == True)).T
    
    X = PtClouds[:,0]
    Y = PtClouds[:,1]
    Z = PtClouds[:,2]
    
    PtCloudPts = np.size(PtClouds,0)
    
    vert_ObjPtCloud = np.zeros((8*PtCloudPts,3));
    
    #get 8 vertices of each pixel
    for q in range(0,PtCloudPts):
        vert_ObjPtCloud[(8*q)+0, 0] = X[q]
        vert_ObjPtCloud[(8*q)+1, 0] = X[q]+1
        vert_ObjPtCloud[(8*q)+2, 0] = X[q]
        vert_ObjPtCloud[(8*q)+3, 0] = X[q]+1
        vert_ObjPtCloud[(8*q)+4, 0] = X[q]
        vert_ObjPtCloud[(8*q)+5, 0] = X[q]+1
        vert_ObjPtCloud[(8*q)+6, 0] = X[q]
        vert_ObjPtCloud[(8*q)+7, 0] = X[q]+1
        
        vert_ObjPtCloud[(8*q)+0, 1] = Y[q]
        vert_ObjPtCloud[(8*q)+1, 1] = Y[q]
        vert_ObjPtCloud[(8*q)+2, 1] = Y[q]+1
        vert_ObjPtCloud[(8*q)+3, 1] = Y[q]+1
        vert_ObjPtCloud[(8*q)+4, 1] = Y[q]
        vert_ObjPtCloud[(8*q)+5, 1] = Y[q]
        vert_ObjPtCloud[(8*q)+6, 1] = Y[q]+1
        vert_ObjPtCloud[(8*q)+7, 1] = Y[q]+1
        
        vert_ObjPtCloud[(8*q)+0, 2] = Z[q]
        vert_ObjPtCloud[(8*q)+1, 2] = Z[q]
        vert_ObjPtCloud[(8*q)+2, 2] = Z[q]
        vert_ObjPtCloud[(8*q)+3, 2] = Z[q]
        vert_ObjPtCloud[(8*q)+4, 2] = Z[q]+1
        vert_ObjPtCloud[(8*q)+5, 2] = Z[q]+1
        vert_ObjPtCloud[(8*q)+6, 2] = Z[q]+1
        vert_ObjPtCloud[(8*q)+7, 2] = Z[q]+1
        
    # number of elements 
    NElmnt = np.sum(V)

    Vertices = vert_ObjPtCloud    
    Vert_rows = np.size(Vertices,0)
    
    # Initialize Nodes Matrix & populate with vertices
    Nodes = np.zeros((Vert_rows,4), dtype=int)
    Indices = np.arange(0,Vert_rows)
    Nodes[:,0] = Indices
    Nodes[:,1:4] = Vertices
    
    
    [NodesMapped, UniqueNodesMapped] = unique(Nodes)
    
    UniqueElements = np.zeros((NElmnt, 9))
    
    NodesMapped = np.vsplit(NodesMapped,NElmnt)
    NodesMapped = np.asarray(NodesMapped, dtype=int)
    
    structure = np.array([0,1,3,2,4,5,7,6])
    NodesMapped = NodesMapped[:,structure,:]
    
    UniqueElements[:,0] = np.arange(1,np.size(UniqueElements,axis=0)+1)
    UniqueElements[:,1:10] = NodesMapped[:,:,0]+1
    
    UniqueNodesMapped[:,1:4] = np.multiply(spacing,UniqueNodesMapped[:,1:4])
    UniqueNodesMapped[:,1] = UniqueNodesMapped[:,1]-np.amin(UniqueNodesMapped[:,1])
    UniqueNodesMapped[:,2] = UniqueNodesMapped[:,2]-np.amin(UniqueNodesMapped[:,2])
    UniqueNodesMapped[:,3] = UniqueNodesMapped[:,3]-np.amin(UniqueNodesMapped[:,3])
    
    return UniqueNodesMapped, UniqueElements

def nodesSet(UniqueNodesMapped, location, strain):
    
    if location == 'EastWest':
        axis = 1
    elif location == 'NorthSouth':
        axis = 2
    else:
        axis = 3
        
    baseCoord = np.min(UniqueNodesMapped[:,axis])    
    baseNodes = np.where((UniqueNodesMapped[:,axis] == baseCoord))    
    baseNodes = np.add(baseNodes,1)

    formatLength = np.divmod(np.size(baseNodes),16)
    
    baseNodesMain = baseNodes[0,0:formatLength[0]*16]
    baseNodesMain = np.reshape(baseNodesMain,(formatLength[0],16))
    baseNodesFinal = baseNodes[0,formatLength[0]*16:]
    
    dispCoord = np.max(UniqueNodesMapped[:,axis])
    
    dispNodes = np.where((UniqueNodesMapped[:,axis] == dispCoord))
    
    dispNodes = np.add(dispNodes,1)
    
    formatLength = np.divmod(np.size(dispNodes),16)
    
    dispNodesMain = dispNodes[0,0:formatLength[0]*16]
    dispNodesMain = np.reshape(dispNodesMain,(formatLength[0],16))
    dispNodesFinal = dispNodes[0,formatLength[0]*16:]
    
    # STRAIN/STRAIN RATE
    
    # axis length
    L=dispCoord-baseCoord
    # applied displacement
    dL = np.round(strain*L,5)      
    
    return baseNodesMain, baseNodesFinal, dispNodesMain, dispNodesFinal, dL


def writeVTK(outPath , caseID, UniqueNodesMapped, UniqueElements):
    
    
    file = path.join(outPath , caseID.replace(".mhd", ".vtk"))
    
    UniqueElements[:,0] = 8
    UniqueElements[:,1:10] = UniqueElements[:,1:10] - 1
    cell_type = np.ones(np.size(UniqueElements,0))*12
    
    with open(file, 'w') as outfile:
        outfile.write('# vtk DataFile Version 3.0\n')
        outfile.write('%s'%date.today()+'\n')
        outfile.write('ASCII\n')
        outfile.write('\n')
        outfile.write('DATASET UNSTRUCTURED_GRID\n')
        outfile.write('POINTS '+'%i'%np.size(UniqueNodesMapped,0)+' double\n')
        np.savetxt(outfile,UniqueNodesMapped[:,1:4],fmt='%11.8e',
                   delimiter=' ', newline='\n')
        outfile.write('\n')
        outfile.write('CELLS '+'%i'%np.size(UniqueElements,0)
                      +' %i'%(9*np.size(UniqueElements,0))+'\n')
        np.savetxt(outfile,UniqueElements,fmt='%i',
                   delimiter=' ', newline='\n')
        outfile.write('\n')
        outfile.write('CELL_TYPES '+'%i'%np.size(UniqueElements,0)+'\n')
        np.savetxt(outfile,cell_type,fmt='%i',
                   delimiter='\n', newline='\n')
        outfile.write('\n')

    
    outfile.close()

if __name__ == '__main__':
    
    #  Read arguments
    inp, write, binary, im, tensor = getInputs()
    
    case = glob(path.join(inp) + im)[0]
    pdir = path.dirname(case)
    caseID = path.split(case)[-1]
    res = ''.join(sorted((set(caseID) ^ set(im)), key = caseID.index))
    outPath = path.join(write)
  
    # Create a new directory if it does not exist 
    isExist = os.path.exists(outPath)
    if not isExist:
        os.makedirs(outPath)
        print("New output directory is created!")
    
               
    # Read image  
    if binary == True:
        [V,spacing]= mhdReader(path.join(inp, caseID))      
    else:
        [V, I,spacing] = threshold(path.join(inp, caseID))
        
   # determine BV/TV  
    bvtv = compute_bvtv(V)
    
    # MIL or GST as fabric tensor M
    if tensor == 'mil':
        # determine mean intercept length tensor M
        [M, info] = compute_mil_tensor(V, spacing, 1000)    
    elif tensor == 'gst':
        # determine mean intercept length tensor M
        [M, info] = compute_gst_tensor(I, spacing, sigma_e=0.1, sigma_i=0.3)
    else:
        print("No known fabric tensor option selected. I am exciting!")
        exit(0)
    
    # writing images for illustration
    # writing tensor as VTK file
    mesh = tensor_to_vtk(M, path.join(outPath, (caseID + "_tensor.vtm")))

    """
    # Quick plot in Python (optional)
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, color="orange", opacity=0.8, show_edges=False)
    plotter.show()
    """ 
    
    # Convert image to hexahedral mesh
    [UniqueNodesMapped, UniqueElements] = cubicMesh(V,np.round(spacing, 5))
    
    # Write mesh to vtk file
    writeVTK(outPath, caseID, UniqueNodesMapped, UniqueElements)
    