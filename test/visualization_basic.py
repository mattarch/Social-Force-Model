import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider

import sys
filename = sys.argv[-1]
f = open(filename, 'rb')
fileAsList = f.readlines()

num_constants = int(fileAsList[0].split()[0])
offset_c = num_constants + 1

AVG_SPEED           = float(fileAsList[1].split()[1])
MAX_SPEED           = float(fileAsList[2].split()[1])
RELAX_TIME          = float(fileAsList[3].split()[1])
WALK_WAY_LENGTH     = float(fileAsList[4].split()[1])
WALK_WAY_WIDTH      = float(fileAsList[5].split()[1])
NUMBER_OF_PEOPLE    = int(fileAsList[6].split()[1])
N_BORDERS           = int(fileAsList[7].split()[1])
TIMESTEP            = float(fileAsList[8].split()[1])
N_TIMESTEP          = int(fileAsList[9].split()[1])
V_ALPHA_BETA        = float(fileAsList[10].split()[1])
SIGMA               = float(fileAsList[11].split()[1])
U_ALPHA_B           = float(fileAsList[12].split()[1])
R                   = float(fileAsList[13].split()[1])
DELTA_T             = float(fileAsList[14].split()[1])
PSI                 = float(fileAsList[15].split()[1])
INFLUENCE           = float(fileAsList[16].split()[1])

fig, axs = plt.subplots(1, 2)

#-------------------------------left plot--------------------------------
# add constants to the plot
textstr = '\n'.join((
    "AVG_SPEED %.2f" % (AVG_SPEED, ),
    "MAX_SPEED %.2f" % (MAX_SPEED, ),
    "RELAX_TIME %.2f" % (RELAX_TIME, ),
    "WALK_WAY_LENGTH %.2f" % (WALK_WAY_LENGTH, ),
    "WALK_WAY_WIDTH %.2f" % (WALK_WAY_WIDTH, ),
    "NUMBER_OF_PEOPLE %.2f" % (NUMBER_OF_PEOPLE, ),
    "N_BORDERS %.2f" % (N_BORDERS, ),
    "TIMESTEP %.2f" % (TIMESTEP, ),
    "N_TIMESTEP %.2f" % (N_TIMESTEP, ),
    "V_ALPHA_BETA %.2f" % (V_ALPHA_BETA, ),
    "SIGMA %.2f" % (SIGMA, ),
    "R %.2f" % (R, ),
    "DELTA_T %.2f" % (DELTA_T, ),
    "PSI %.2f" % (PSI, ),
    "INFLUENCE %.2f" % (INFLUENCE, ),
))

axs[0].text(0.5, 0.5, textstr, transform=axs[0].transAxes, fontsize=12,verticalalignment='center',horizontalalignment='center')
axs[0].axis('off')
axs[1].axis('equal')

categories = []

color1=(0.69411766529083252, 0.3490196168422699, 0.15686275064945221, 1.0)
color2=(0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)

colormap = np.array([color1,color2])

#-------------------------------right plot--------------------------------

# plot borders for simple scenario
axs[1].plot([0,WALK_WAY_LENGTH],[WALK_WAY_WIDTH,WALK_WAY_WIDTH],'k')
axs[1].plot([0,WALK_WAY_LENGTH],[0,0],'k')

# plot inital values
x_initial = []
y_initial = []
u_initial = []
v_initial = []
x_dest = []
y_dest = []
for i in range(NUMBER_OF_PEOPLE):
    if(i % 2):
        categories.append(0)
    else:
        categories.append(1)

    line = fileAsList[i+offset].split()
    speed = 1 #float(line[2])
    x_initial.append(float(line[0])) #x position
    y_initial.append(float(line[1])) #y position
    u_initial.append(speed * float(line[3])) 
    v_initial.append(speed * float(line[4]))
    x_dest.append(float(line[5]))
    y_dest.append(float(line[6]))

# plot inital points
scat = axs[1].scatter(x_initial,y_initial,marker='o', c=colormap[categories])

# plot destination points
scat_dest = axs[1].scatter(x_dest,y_dest,marker='o', color="blue")

# plot direct_x, direct_y arrows
quiv = axs[1].quiver(x_initial,y_initial,u_initial,v_initial)


# Slider
axamp = plt.axes([0.6, .03, 0.25, 0.02])

samp = Slider(axamp, 'Timesteps', 0, N_TIMESTEP, valinit=0, valstep=1,)

#-------------------------------Animation--------------------------------

# Animation controls
is_manual = True            # True if user has taken control of the animation
interval = 500              # ms, time between animation frames
loop_len = 0.5 * N_TIMESTEP # seconds per loop

def read_data(iter):
    offset = int(NUMBER_OF_PEOPLE * iter) + offset_c
    x = []
    y = []
    u = []
    v = []
    speed = 1.0
    for i in range(NUMBER_OF_PEOPLE):
        line = fileAsList[i+offset].split()
        x.append(float(line[0])) #x position
        y.append(float(line[1])) #y position
        u.append(speed * float(line[7])) 
        v.append(speed * float(line[8]))
    return x, y, u, v
   
         


# slider functions   

def update_slider(val):
    global is_manual
    is_manual=True
    update(val)

def update(val):
    global scat
    global quiv
    if val < N_TIMESTEP:
        x, y, u, v = read_data(val)
        scat.set_offsets(np.c_[x, y])
        quiv.set_offsets(np.c_[x, y])
        quiv.set_UVC(u,v)

    # redraw canvas while idle
    fig.canvas.draw_idle()

def update_plot(num):
    global is_manual
    if is_manual:
        return fig, # don't change

    val = (samp.val + 1) % samp.valmax
    samp.set_val(val)
    is_manual = False
    return fig,

def on_click(event):
    # Check where the click happened
    (xm,ym),(xM,yM) = samp.label.clipbox.get_points()
    if xm < event.x < xM and ym < event.y < yM:
        # Event happened within the slider, ignore since it is handled in update_slider
        return
    else:
        # user clicked somewhere else on canvas = unpause
        global is_manual
        is_manual=False

# call update function on slider value change
samp.on_changed(update_slider)

fig.canvas.mpl_connect('button_press_event', on_click)

ani = animation.FuncAnimation(fig, update_plot, interval=interval,frames=100)

plt.show()
