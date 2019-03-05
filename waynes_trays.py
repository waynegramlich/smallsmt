#!/usr/bin/env python

#<----------------------------------------- 100 characters --------------------------------------->|

from small_smt import *

def main():
    # Create *ezcad* 3.0:
    ezcad = EZCAD3(0)

    # Define all of the channels:
    full_length   = True
    normal_length = False
    edge_depth = L(mm=0.6)
    channel8f  = Channel(8.0,  full_length,   L(mm=0.64),  4.55,  2.40, 1.75 +  3.50)
    channel8n  = Channel(8.0,  normal_length, L(mm=0.64),  4.55,  2.40, 1.75 +  3.50)
    channel12n = Channel(12.0, normal_length, edge_depth,  8.20,  6.40, 1.75 +  5.50)
    channel12f = Channel(12.0, full_length,   edge_depth,  8.20,  6.40, 1.75 +  5.50)
    channel16  = Channel(16.0, normal_length, edge_depth, 12.10,  7.90, 1.75 +  7.50)
    channel24  = Channel(24.0, normal_length, edge_depth, 20.10, 11.90, 1.75 + 11.50)

    # Channels are named using the last digit of the channel width:
    #        8mm  => 8
    #        12mm => 2
    #        16mm => 6
    #        20mm => 0
    #        24mm => 4

    # 1 12mm channel => "2":
    colors_2   = ("aqua", "medium_purple", "purple")
    channels_2 = (channel12n, )
    tray_2     = Tray_Specification("Tray_2", 0, colors_2, channels_2)

    # 5 8mm channels => "88888" :
    colors_88888   = ("light_blue", "blue", "deep_sky_blue")
    channels_88888 = (channel8n, channel8n, channel8f, channel8n, channel8n)
    tray_88888     = Tray_Specification("Tray_88888", 1, colors_88888, channels_88888)

    # 12mm/8mm/12mm channels => "282":
    colors_282   = ("pink", "red", "dark_red")
    channels_282 = (channel12n, channel8f, channel12n)
    tray_282     = Tray_Specification("Tray_282", 1, colors_282, channels_282)

    # Create *tray_specifications* which is a tuple of the *Tray_Specifcation* objects:
    tray_specifications = (tray_2, )
    #tray_specifications = (tray_2, tray_88888, tray_282)

    # Create the top level *small_smt* object using *tray_specifations*.
    # To disable the visualization cut-outs, set `debug` to `False`:
    small_smt = SmallSMT(None, "SmallSMT", tray_specifications, debug=True)

    # Cause the system to be generated:
    small_smt.process(ezcad)

if __name__ == "__main__":
    main()
