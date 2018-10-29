import sys
import os
import datetime
import imageio


def make_gif(png_dir, name, loops, fps):
    '''
    INPUTS:
    loop : int
    The number of iterations. Default 0 (meaning loop indefinitely).

    duration : {float, list}
    The duration (in seconds) of each frame.
    Either specify one value that is used for all frames, or one value for
    each frame. Note that in the GIF format the duration/delay is expressed in
    hundredths of a second, which limits the precision of the duration.

    fps : float
    The number of frames per second. If duration is not given, the duration
    for each frame is set to 1/fps. Default 10.

    palettesize : int
    The number of colors to quantize the image to. Is rounded to the nearest
    power of two. Default 256.

    subrectangles : bool
    If True, will try and optimize the GIF by storing only the rectangular
    parts of each frame that change with respect to the previous.
    Default False.
    '''

    images = []
    temp = []
    for file_name in os.listdir(png_dir):
        if file_name.endswith('.png'):
            temp.append(file_name)

    #temp = sorted(temp, key=lambda x: x[2], reverse=True)
    # print(temp)

    for file_name in temp:
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
        print(file_name)

    output_file = png_dir + name
    imageio.mimsave(output_file, images, loop=loops, fps=fps)

directory = '/Users/arpj/Dropbox/Documentos/ETH Zuerich/ETH MAS Thesis/presentations/presentation3/diagrams/rc/b/'
name = 'rc.gif'
make_gif(directory, name, 1, 24)
