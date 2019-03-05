# SmallSMT Pick-and-Place Machine

Parts designed for a SmallSMT pick-and-place machine.

## Introduction

The company [SmallSMT](https://www.smallsmt.biz/) sells a desktop pick and place
machine for small production runs of PCB's (Printed Circuit Boards.)
This reposisitory has designs for some cut tape part holding trays for their
mid-range VP2-2000S/WE, VP-2500DP/WN, and VP-2500DP/WNE machines.

These designs use some custom Python code for designing parts called EZCAD3.
[EZCAD3](https://github.com/waynegramlich/ezcad3) is basically Python code that
layers on top of the `openscad` program to design parts.  In addition, it has
CNC tool path generation.

## Installation

In order to use this code please do the following:

        cd {somewhere}
	git clone https://github.com/waynegramlich/smallsmt.git
	git clone https://github.com/waynegramlich/ezcad3.git
	sudo -H apt-get install openscad
	sudo -H apt-get install view3dscene

## Execution

To generate the parts, do the following:

        cd {somewhere}/smallsmt
	./waynes_trays.py   # This takes a while!
	view3dscene wrl/SmallSMT.wrl  # Use to view the generated stuff.

Please copy `waynes_trays.py` to another file and edit it to suit your needs.

