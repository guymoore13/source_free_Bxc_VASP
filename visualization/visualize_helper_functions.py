 
import numpy as np
import numpy.linalg as npla

import pandas as pd

from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import VolumetricData

import plotly
import plotly.graph_objects as go

def coarsen_voldata(phi, coarseness):
    # From: https://stackoverflow.com/questions/34689519/how-to-coarser-the-2-d-array-data-resolution
    temp = phi.reshape((
        phi.shape[0] // coarseness, coarseness,
        phi.shape[1] // coarseness, coarseness,
        phi.shape[2] // coarseness, coarseness,
    ))
    coarse_phi = np.mean(temp, axis=(1,3,5))
    
    return coarse_phi

def get_lat_mat(struct):
    
    lat_mat = np.array(struct.lattice.matrix.T)
    # lat_mat = np.array(struct.lattice.matrix)
    # lat_mat = np.eye(3)
    
    # lat_mat /= npla.det(lat_mat)**(1/3)
    
    return lat_mat

def get_plot_data_vector(
    voldata, struct, 
    coarseness=2,
):
    
    lat_mat = get_lat_mat(struct)
    
    diff_x_coarse = coarsen_voldata(voldata.data['diff_x'], coarseness)
    diff_y_coarse = coarsen_voldata(voldata.data['diff_y'], coarseness)
    diff_z_coarse = coarsen_voldata(voldata.data['diff_z'], coarseness)
    
    nshape = diff_x_coarse.shape
    nnn = np.prod(nshape)
    
    xf,yf,zf = np.linspace(0.0,1.0,nshape[0]), np.linspace(0.0,1.0,nshape[1]), np.linspace(0.0,1.0,nshape[2])
    xxf,yyf,zzf = np.meshgrid(xf,yf,zf, indexing='ij')
    
    xyzf = np.vstack([np.reshape(xxf,[nnn]), np.reshape(yyf,[nnn]), np.reshape(zzf,[nnn])])
    xyz = lat_mat @ xyzf
    xx,yy,zz = xyz[0,:], xyz[1,:], xyz[2,:]
    
    # xx,yy,zz = xx.reshape([nnn]), yy.reshape([nnn]), zz.reshape([nnn])
    diff_x_coarse = diff_x_coarse.reshape([nnn])
    diff_y_coarse = diff_y_coarse.reshape([nnn])
    diff_z_coarse = diff_z_coarse.reshape([nnn])
    
    ux, uy, uz = voldata.data['diff_x'], voldata.data['diff_y'], voldata.data['diff_z']
    # ux, uy, uz = diff_x_coarse, diff_y_coarse, diff_z_coarse
    norm = np.sqrt(ux**2 + uy**2 + uz**2)
    # norm = np.mean(norm)
    norm = np.max(norm)
    
    print("Dimensions:", nshape, np.min(nshape))
    print("Norm = ", norm)
    
    print("Integral (of sorts) x,y,z components = %f %f %f" % (np.mean(ux), np.mean(uy), np.mean(uz)))
    
    return nshape, xx, yy, zz, diff_x_coarse, diff_y_coarse, diff_z_coarse

def get_plot_data_stream(
    nshape, struct, 
    stream_coarseness=4,
):
    
    lat_mat = get_lat_mat(struct)
    
    snshape = [n//stream_coarseness for n in nshape]
    snnn = np.prod(snshape)
    
    xf,yf,zf = np.linspace(0.0,1.0,snshape[0]), np.linspace(0.0,1.0,snshape[1]), np.linspace(0.0,1.0,snshape[2])
    xxf,yyf,zzf = np.meshgrid(xf,yf,zf, indexing='ij')
    
    xyzf = np.vstack([np.reshape(xxf,[snnn]), np.reshape(yyf,[snnn]), np.reshape(zzf,[snnn])])
    xyz = lat_mat @ xyzf
    xxs,yys,zzs = xyz[0,:], xyz[1,:], xyz[2,:]
        
    snnn = np.prod(snshape)
    xxs,yys,zzs = xxs.reshape([snnn]), yys.reshape([snnn]), zzs.reshape([snnn])
    
    print("Stream dimensions:", snshape, np.min(snshape))
    
    return xxs, yys, zzs
    

def visualize_voldata(
    filename_a, filename_b, 
    coarseness, stream_coarseness,
    compute_torque=False):
    
    # read in vector field data
    poscar_a, data_a, data_aug_a = VolumetricData.parse_file(filename_a)
    struct_a = poscar_a.structure
    voldata_a = VolumetricData(structure=struct_a, data=data_a)
    
    nshape_a, xx_a, yy_a, zz_a, diff_x_coarse_a, diff_y_coarse_a, diff_z_coarse_a = \
        get_plot_data_vector(voldata_a, struct_a, coarseness)
    
    # read in stream line data
    poscar_b, data_b, data_aug_b = VolumetricData.parse_file(filename_b)
    struct_b = poscar_b.structure
    voldata_b = VolumetricData(structure=struct_b, data=data_b)
    
    nshape_b, xx_b, yy_b, zz_b, diff_x_coarse_b, diff_y_coarse_b, diff_z_coarse_b = \
        get_plot_data_vector(voldata_b, struct_b, coarseness)
    xxs_b, yys_b, zzs_b = get_plot_data_stream(nshape_b, struct_b, stream_coarseness)
    
    # compute torque
    ashape = voldata_a.data['diff_x'].shape
    tau_xc = {d:np.zeros(ashape) for d in ['total', 'diff_x', 'diff_y', 'diff_z']}
    
    if compute_torque:
        for i in range(ashape[0]):
            for j in range(ashape[1]):
                for k in range(ashape[2]):
                    av = np.array([voldata_a.data['diff_x'][i,j,k], voldata_a.data['diff_y'][i,j,k], voldata_a.data['diff_z'][i,j,k]])
                    bv = np.array([voldata_b.data['diff_x'][i,j,k], voldata_b.data['diff_y'][i,j,k], voldata_b.data['diff_z'][i,j,k]])
                    tv = np.cross(av, bv)
                    tau_xc['diff_x'][i,j,k], tau_xc['diff_y'][i,j,k], tau_xc['diff_z'][i,j,k] = tv
        voldata_t = VolumetricData(structure=struct_a, data=tau_xc)
        nshape_t, xx_t, yy_t, zz_t, diff_x_coarse_t, diff_y_coarse_t, diff_z_coarse_t = \
            get_plot_data_vector(voldata_t, struct_a, coarseness)
    
    struct = struct_a.copy()
    
    lat_mat = get_lat_mat(struct)
    
    data = []
    
    data.append(go.Cone(
        x=xx_a, y=yy_a, z=zz_a,
        u=diff_x_coarse_a,
        v=diff_y_coarse_a,
        w=diff_z_coarse_a,
        colorscale='Inferno',
        showscale = False,
        sizeref=1.5,
    ))
    
    data.append(go.Streamtube(
        x=xx_b, y=yy_b, z=zz_b,
        u=diff_x_coarse_b,
        v=diff_y_coarse_b,
        w=diff_z_coarse_b,
        starts = dict(
            x = xxs_b, y = yys_b, z = zzs_b,
        ),
        sizeref = 0.2,
        colorscale = 'Inferno',
        showscale = False,
        maxdisplayed = 1000,
        opacity = 0.1,
    ))
    
    # create box
    a, b, c = struct.lattice.a, struct.lattice.b, struct.lattice.c
    xxf = [[0,1],[0,0],[0,0],[0,1],[1,1],[1,1], [0,0],[0,0],[1,0],[1,1],[1,0],[1,1]]
    yyf = [[0,0],[0,1],[0,0],[1,1],[0,1],[1,1], [1,0],[1,1],[0,0],[0,0],[1,1],[1,0]]
    zzf = [[0,0],[0,0],[0,1],[1,1],[1,1],[0,1], [1,1],[1,0],[1,1],[1,0],[0,0],[0,0]]
    for xf,yf,zf in zip(xxf,yyf,zzf):
        xyz = lat_mat @ np.vstack([xf,yf,zf])
        data.append(go.Scatter3d(
            x=xyz[0,:], y=xyz[1,:], z=xyz[2,:],
            mode='lines',
            marker=None,
            line=dict(
                color='grey',
                width=2
            )
        ))
    
    fig = go.Figure(data=data)
    
    fig.update_layout(
        width=800, height=800,
        autosize=False,
        margin = {'l':0,'r':0,'t':0,'b':0},
        scene=dict(
            aspectratio=dict(x=1, y=1, z=1),
            # camera_eye=dict(x=1.0, y=1.0, z=1.0),
            camera_eye=dict(
                x=0.333*struct_a.lattice.a, 
                y=0.333*struct_a.lattice.b, 
                z=0.333*struct_a.lattice.c),
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        font=dict(
            family="Serif",
            size=24,
        ),
        showlegend=False,
    )
        
    return fig, tau_xc

