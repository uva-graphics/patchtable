
import sys
import skimage
import skimage.io

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'greyscale in.png out.png'
        print >> sys.stderr, 'Converts image to greyscale, output is 3-channel RGB'
        sys.exit(1)
        
    a = skimage.img_as_float(skimage.io.imread(args[0]))
    L = a[:,:,0]*0.30 + a[:,:,1]*0.59 + a[:,:,2]*0.11
    for channel in range(3):
        a[:,:,channel] = L
    skimage.io.imsave(args[1], a)
    
if __name__ == '__main__':
    main()

