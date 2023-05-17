import pandas as pd
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
from datetime import datetime
import time
from PIL import Image

def create_image(Data, Y, X):
    raw_image = np.full((Y,X,4), 255)
    H = np.max(Data)
    L = np.min(Data)
    
    raw_image[:,:,1] = (Data[:,:]-L)*255/(H-L)
    raw_image[:,:,0] = 0
    raw_image[:,:,2] = 0
    raw_image[:,:,3] = 255
    return raw_image

start_time = datetime.now()

H = 1
l = 0.3
L = 0.7
X = 1
T = 1
Kurant = 0.1

def calc_sigma (ts, xs, u):
	return ts/xs*np.max(u)

def calc_initial_t (t):
	if t < l or t >= L:
		return 0
	if t >= l and t < (l+L)/2:
		return (2*H/(L-l))*(t-l)
	if t >= (l+L)/2 and t < L:
		return -(2*H/(L-l))*(t-L)

def calc_initial_x (x):
	if x < l or x >= L:
		return 0
	if x >= l and x < L:
		return H

def left_angel(ud, udl, Kurant):
	uds = ud**2
	udls = udl**2
	return Kurant * (uds - udls)/(2) + ud

def scheme_CIR(ud, udl, udr, Kurant):
    if ud < 0:
        return ud - (udr**2 - ud**2)*Kurant/(2)
    else:
        return ud - (ud**2 - udl**2)*Kurant/(2 )

def create_GIF(ts, Kurant):

	dir_name = "./temp/"
	os.mkdir(dir_name)

	xs = ts * int(1 / Kurant)

	Nx = int(X/xs)+1
	Nframes = int(T/ts)+1
 
	if (Nframes > 100):
		stepes = int(Nframes / 100)
	else:
		stepes = 1

	# CALCULATION
	u = np.full((Nframes,Nx), 0, dtype="float64")
	u[0,:] = np.linspace(0,X,Nx)
	u[:,0] = np.linspace(0,T,Nframes)
	for i in range (0, Nx):
		u[0,i] = calc_initial_x(u[0,i])
	for i in range (0, Nframes):
		u[i,0] = calc_initial_t(u[i,0])

	CIR = u.copy()

	for i in range (1, Nframes):
		for j in range (1,Nx):
			if (j != Nx-1):
				CIR[i,j] = scheme_CIR(CIR[i-1,j],CIR[i-1,j-1],CIR[i-1,j+1],Kurant)
			else:
				CIR[i,j] = scheme_CIR(CIR[i-1,j],CIR[i-1,j-1], CIR[i-1, j],Kurant)

	# MAKE FRAMES
	x = np.linspace(0, X, Nx)

	for frame in range(0, Nframes, stepes):
		fig = plt.figure(figsize = [12,7])
		fig = plt.plot(x, CIR[frame, :], '-b', ms = 0.1)
		fig = plt.xlabel('x')
		fig = plt.xlim([0, X])
		fig = plt.ylim([0, 1.2])
		fig = plt.ylabel('u(x)')
		fig = plt.title(f"t = {round(frame * ts, 4)}")
		fig = plt.savefig(dir_name + f"{frame}.png")
		plt.close(fig)

	# GIF CREATION
	im = []
	for frame in range(0, Nframes, stepes):
		img = Image.open(dir_name + f"{frame}.png")
		img.load()
		im.append(img) 

	im[0].save(f"./GIF/CIR_{xs}.gif", save_all = "True", append_images = im[1:], duration = 10)

def do_CIR(ts, Kurant):

	xs = ts * int(1 / Kurant)
	print("CIR prec: ", xs)

	Nx = int(X/xs)+1
	Nt = int(T/ts)+1

	u = np.full((Nt,Nx), 0, dtype="float128")
	u[0,:] = np.linspace(0,X,Nx)
	u[:,0] = np.linspace(0,T,Nt)
	for i in range (0, Nx):
		u[0,i] = calc_initial_x(u[0,i])
	for i in range (0, Nt):
		u[i,0] = calc_initial_t(u[i,0])

	CIR = u.copy()

	for i in range (1, Nt):
		for j in range (1,Nx):
			if (j != Nx-1):
				CIR[i,j] = scheme_CIR(CIR[i-1,j],CIR[i-1,j-1],CIR[i-1,j+1],Kurant)
			else:
				CIR[i,j] = scheme_CIR(CIR[i-1,j],CIR[i-1,j-1], CIR[i-1, j],Kurant)
	# print(CIR)

	calc_time = datetime.now() - start_time
	print("calc: ", calc_time)

	df = pd.DataFrame(data=CIR)
	df.to_csv(f"./tables/CIR-{xs}.csv")

	save_time = datetime.now() - start_time - calc_time
	print("save: ", save_time)

	x = np.linspace(0, X, Nx)
	t = np.linspace(0, T, Nt)

	raw = create_image(CIR,Nt,Nx)
	im = Image.fromarray(raw.astype(np.uint8))
	size = (1000,1000)
	im = im.resize(size)
	im = im.transpose(Image.Transpose.FLIP_TOP_BOTTOM)
	im.save(f"./images/CIR-{xs}.png")

	total_time = datetime.now() - start_time - save_time - calc_time
	print("show: ", total_time)

def do_CIR(ts, Kurant):

    xs = ts * int(1 / Kurant)
    print("CIR prec: ", xs)

    Nx = int(X/xs)+1
    Nt = int(T/ts)+1

    u = np.full((Nt,Nx), 0, dtype="float128")
    u[0,:] = np.linspace(0,X,Nx)
    u[:,0] = np.linspace(0,T,Nt)
    for i in range (0, Nx):
        u[0,i] = calc_initial_x(u[0,i])
    for i in range (0, Nt):
        u[i,0] = calc_initial_t(u[i,0])

    CIR = u.copy()

    for i in range (1, Nt):
        for j in range (1,Nx):
            if (j != Nx-1):
                CIR[i,j] = scheme_CIR(CIR[i-1,j],CIR[i-1,j-1],CIR[i-1,j+1],Kurant)
            else:
                CIR[i,j] = scheme_CIR(CIR[i-1,j],CIR[i-1,j-1], CIR[i-1, j],Kurant)
    # print(CIR)

    calc_time = datetime.now() - start_time
    print("calc: ", calc_time)

    df = pd.DataFrame(data=CIR)
    df.to_csv(f"./tables/CIR-{xs}.csv")

    save_time = datetime.now() - start_time - calc_time
    print("save: ", save_time)

    x = np.linspace(0, X, Nx)
    t = np.linspace(0, T, Nt)

    raw = create_image(CIR,Nt,Nx)
    im = Image.fromarray(raw.astype(np.uint8))
    size = (1000,1000)
    im = im.resize(size)
    im = im.transpose(Image.Transpose.FLIP_TOP_BOTTOM)
    im.save(f"./images/CIR-{xs}.png")

    total_time = datetime.now() - start_time - save_time - calc_time
    print("show: ", total_time)

t = 0.0005
Kurant = 0.1
do_CIR(t, Kurant)