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
    channel12f = Channel(12.0, full_length,   edge_depth,  8.20,  6.40, 1.75 +  5.50)
    channel12n = Channel(12.0, normal_length, edge_depth,  8.20,  6.40, 1.75 +  5.50)
    channel16n = Channel(16.0, normal_length, edge_depth, 12.10,  7.90, 1.75 +  7.50)
    channel20n = Channel(20.0, normal_length, edge_depth, 12.10,  7.90, 1.75 +  7.50)
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

    # 1 16mm channel => "6":
    colors_6   = ("blanched_almond", "brown", "chocolate")
    channels_6 = (channel20n, )
    tray_6     = Tray_Specification("Tray_6", 0, colors_6, channels_6)

    # 6 8mm channels => "888888" :
    colors_888888   = ("light_blue", "blue", "deep_sky_blue")
    channels_888888 = (channel8n, channel8f, channel8f, channel8f, channel8f, channel8f)
    tray_888888     = Tray_Specification("Tray_888888", 2, colors_888888, channels_888888)

    # 16mm/16mm/8mm channels => "668":
    colors_668  = ("pink", "red", "dark_red")
    channels_668 = (channel16n, channel16n, channel8f)
    tray_668     = Tray_Specification("Tray_668", 1, colors_668, channels_668)

    # Create *tray_specifications* which is a tuple of the *Tray_Specifcation* objects:
    #tray_specifications = (tray_2, )
    #tray_specifications = (tray_2, tray_6)
    tray_specifications = (tray_2, tray_6, tray_888888, tray_668)

    # Create the top level *small_smt* object using *tray_specifations*.
    # To disable the visualization cut-outs, set `debug` to `False`:
    small_smt = SmallSMT(None, "SmallSMT", tray_specifications, debug=True)

    # Cause the system to be generated:
    small_smt.process(ezcad)

if __name__ == "__main__":
    main()
