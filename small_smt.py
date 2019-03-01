#!/usr/bin/env python

#<----------------------------------------- 100 characters --------------------------------------->|

# Classes are listed alphabetically.

from EZCAD3 import *

class Channel:
    """ *Channel*: Represents the information about the cut tape channel.
    """

    def __init__(self,
      tape_width, tape_length, tape_edge_depth, pocket_width, pocket_depth, pocket_offset):
	""" *Channel*:  Initialize the *Channel* object (i.e. *self*).

	The arguments are:
	* *tape_width*: Width of the tape in millimeters.
	* *tape_length*: Length of the tap in millimeters.
	* *tape_edge_depth*: Depth of the tape edge in millimeters.
	* *pocket_width*: Width of component pocket in millimeters in Y direction.
	* *pocket_depth*: Depth of component pocket in millimeters.
	* *pocket_offset*: Component pocket Y offset from top tape edge (holes on top).
        If there is no pocket, *pocket_width*, *pocket_depth*, and *pocket_offset* are 0.
	"""

	# Verify argument types:
	assert isinstance(tape_width, float)
	assert isinstance(tape_length, float)
	assert isinstance(tape_edge_depth, float)
	assert isinstance(pocket_width, float)
	assert isinstance(pocket_depth, float)
	assert isinstance(pocket_offset, float)

	# Stuff values into *channel* object (i.e. *self*):
	channel = self
	channel.tape_width      = L(mm=tape_width)
	channel.tape_length     = L(mm=tape_length)
	channel.tape_edge_depth = L(mm=tape_edge_depth)
	channel.pocket_width    = L(mm=pocket_width)
	channel.pocket_depth    = L(mm=pocket_depth)
	channel.pocket_offset   = L(mm=pocket_offset)

class Frame(Part):
    """ *Frame*: Represents the relevant *Frame* portion of pick and place machine.
    """

    def __init__(self, up, name):
	""" *Frame* Initialize the *Frame* object (i.e. *self*):
	"""

	# Standard initialization sequence:
	frame = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(frame, up, name)

	# Create the various portions of the *frame*:
	frame.east_rail_ = Rail(frame, "East_Rail", True)
	frame.west_rail_ = Rail(frame, "West_Rail", False)

    def construct(self):
	""" *Frame*: Construct the *Frame* assembly object (i.e. *self*.)
	"""

	# Specify various frame constants here:
	frame = self
	frame.rail_dx_l        = rail_dx        = L(mm=9.20)
	frame.rail_dx_pitch_l  = rail_dx_pitch  = L(mm=216.0) + rail_dx
	frame.rail_dy_l        = rail_dy        = L(mm=325.0)
	frame.holes_pitch_dy_l = holes_pitch_dy = L(mm=20.00)
	frame.rail_top_z_l     = rail_top_z     = L(mm=-4.00)
	frame.holes_dy_l       = holes_dy     = float(int(rail_dy/holes_pitch_dy)) * holes_pitch_dy

class Rail(Part):
    """ *Rail*: Represents one of the rails.
    """

    def __init__(self, up, name, is_east):
	""" *Rail*: Initialize the *Rail* object (i.e. *self*.)
	"""

	# Standard initialization sequence:
	rail = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(is_east, bool)
	Part.__init__(rail, up, name)

	# Hang onto the *is_east*:
	rail.is_east_b = is_east

    def construct(self):
	""" *Rail*: Construct the *Rail* object (i.e. *self*.)
	"""

	# Grab some values from *frame*:
	rail          = self
	is_east       = rail.is_east_b
	frame         = rail.up
	rail_dx       = frame.rail_dx_l
	rail_dx_pitch = frame.rail_dx_pitch_l
	rail_dy       = frame.rail_dy_l
	rail_top_z    = frame.rail_top_z_l
	holes_dy      = frame.holes_dy_l

	# Load some constants into *rail*:
	rail.dx_l             = dx             = rail_dx
	rail.dy_l             = dy             = rail_dy
	rail.dz_l             = dz             = L(mm=7.00)
	rail.holes_pitch_dy_l = holes_pitch_dy = L(mm=20.00)

	# Compute some X coordinates:
	zero = L()
	x10 =  rail_dx_pitch/2 + dx/2
	x8  =  rail_dx_pitch/2
	x6  =  rail_dx_pitch/2 - dx/2
	x5  =  zero
	x4  = -rail_dx_pitch/2 + dx/2
	x2  = -rail_dx_pitch/2
	x0  = -rail_dx_pitch/2 - dx/2

	# Compute some Y coordinates:
	y10 =  dy/2
	y8  =  holes_dy/2
	y5  =  zero
	y2  = -holes_dy/2
	y0  = -dy/2

	# Compute some X coordinates:
	z10 = rail_top_z
	z5  = rail_top_z - dz/2
	z0  = rail_top_z - dz

	# Create the Rail as a block:
	material = Material("Plastic", "HDPE")
	color = Color("brown")
	if is_east:
	    x_west = x6
	    x_center = x8
	    x_east = x10
	else:
	    x_west = x0
	    x_center = x2
	    x_east = x4

	# Start with a block of *material*:
	corner1 = P(x_west,  y0,  z0)
	corner2 = P(x_east, y10, z10)
	rail.block("Block", material, color, corner1, corner2, "")
	rail.vice_mount("Top_Vice", "t", "n", "", zero, zero)	

	index = 0
	y = y2
	while y <= y8:
	    start = P(x_center, y, z10)
	    stop1  = P(x_center, y, z5)
	    stop2  = P(x_center, y, z0)
	    comment1 = "Hole {0}a".format(index)
	    comment2 = "Hole {0}2".format(index)
	    diameter1 = dx/2
	    diameter2 = L(inch=0.0890) # #2-56:close
	    rail.hole(comment1, diameter1, start, stop1, "f")
	    rail.hole(comment2, diameter2, start, stop2, "t")
	    y += holes_pitch_dy

	# Stuff some values into *rail*:
	rail.wide_hole_dz_l         = z10 - z5
	rail.wide_hole_diameter_l   = diameter1
	rail.narrow_hole_diameter_l = diameter2

class SmallSMT(Part):
    """ *SmallSMT*: Represents the entir small SmallSMT pick-and-place machine of interest.
    """

    def __init__(self, up, name):
	""" *SmallSMT*: Initialize the *SmallSMT* object (i.e. *self*.)
	"""

	# Standard initialization sequence:
	small_smt = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(small_smt, up, name)

	small_smt.frame_ = Frame(small_smt, "Frame")
	small_smt.trays_ = Trays(small_smt, "Trays")

    def construct(self):
	""" *SmallSMT*: Construct the *SmallSMT* object (i.e. *self8:)
	"""

	pass

class Tray(Part):
    """ *Tray*: Represents one cut-tape tray. """

    def __init__(self, up, name, colors, first_hole, holes_count, channels):
	""" *Tray*: Initialize *Tray* assembly object.
	"""
	
	# Standard initialization sequence:
	tray = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(tray, up, name)

	# Verify additional arguments:
	assert isinstance(colors, tuple) and len(colors) == 3
	assert isinstance(first_hole, int) and first_hole >= 0
	assert isinstance(holes_count, int) and holes_count >= 1
	assert isinstance(channels, tuple) and len(channels) >= 1

	# Define the tray components:
	tray.base_ = Tray_Base(tray, name + "_Base", Color(colors[0]))
	#tray.lip_  = Tray_Lid(tray,  name + "_Lid", Color(colors[1]))
	#tray.top_  = Tray_Top(tray,  name + "_Top", Color(colors[2]))

	# Save the additional arguments into *tray*:
	tray.first_hole_i        = first_hole
	tray.holes_count_i       = holes_count
	tray.channels_o          = channels
	zero                     = L()
	tray.channel_y_centers_o = [zero for channel in channels]

    def construct(self):
	""" *Tray*: Constructs the *Tray* object. """

	# Grab some values from *tray*:
	tray        = self
	first_hole  = tray.first_hole_i
	holes_count = tray.holes_count_i
	channels    = tray.channels_o
	assert isinstance(channels, tuple) and len(channels) > 0

	# Grap some useful *Part*'s:
	trays     = tray.up
	small_smt = trays.up
	frame     = small_smt.frame_
	west_rail = frame.west_rail_
	
	rail_dx        = frame.rail_dx_l
	rail_dx_pitch  = frame.rail_dx_pitch_l
	holes_pitch_dy = frame.holes_pitch_dy_l
	holes_dy       = frame.holes_dy_l

	zero = L()
	tray.shave_dy_l   = shave_dy   = L(inch=0.025)
	tray.dx_l         = dx         = rail_dx_pitch + 2 * rail_dx 
	tray.dy_l         = dy         = (holes_count + 1) * holes_pitch_dy - 2 * shave_dy
	tray.north_y_l    = north_y    = ( -holes_dy/2 +
                                      (first_hole + holes_count + 0.5) * holes_pitch_dy ) - shave_dy
	tray.south_y_l    = south_y   = -holes_dy/2 + (first_hole - 0.5) * holes_pitch_dy + shave_dy
	tray.base_top_z_l = base_top_z = zero   # Top surface of *Tray_Base*

	#print("tape_widths={0}".format(tape_widths))
	tray.channel_widths_o = channel_widths = [channel.tape_width for channel in channels]
	assert isinstance(channel_widths, list)
	zero = L()
	total_channels_width = zero
	for channel_width in channel_widths:
	    #print("channel_width={0:m}".format(channel_width))
	    total_channels_width += channel_width
	#print("total_channels_width={0:m}".format(total_channels_width))
	left_over = dy - total_channels_width
	gap_dy = left_over / (len(channel_widths) + 1)
	#print("dy={0:m} left_over={1:m} gap_dy={2:m}".format(dy, left_over, gap_dy))

	tray.channel_y_centers_o = channel_y_centers = list()
	#print("south_y={0:m}".format(south_y))
	y = south_y + gap_dy
	for channel_width in channel_widths:
	    channel_y_center = y + channel_width/2
	    #print("y={0:m} channel_width={1:m} channel_y_center={2:m}".
	    #  format(y, channel_width, channel_y_center))
	    channel_y_centers.append(channel_y_center)
	    y += channel_width + gap_dy
	#print("north_y={0:m}".format(north_y))

class Tray_Base(Part):
    """ *Tray_Base*: Represents the Tray base that holes the cut tape. """

    def __init__(self, up, name, color):
	""" *Tray_Base*: Represents
	"""

	# Standard initialization sequence:
	tray_base = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(color, Color)
	Part.__init__(tray_base, up, name)
	tray_base.color_o = color

    def construct(self):
	""" *Tray_Base*: Construct the *Tray_Base* object.
	"""

	# Extract some *Part*'s from the *Tray_Base* object (i.e. *self*):
	tray_base = self
	tray      = tray_base.up
	trays     = tray.up
	small_smt = trays.up
	frame     = small_smt.frame_
	west_rail = frame.west_rail_
	east_rail = frame.east_rail_

	# Grab some values from the *Part*'s:
	color                   = tray_base.color_o
	name                    = tray_base.name_get()
	tray_dx                 = tray.dx_l
	tray_dy                 = tray.dy_l
	tray_shave_dy           = tray.shave_dy_l
	tray_north_y            = tray.north_y_l
	tray_south_y            = tray.south_y_l
	tray_base_top_z         = tray.base_top_z_l
	tray_channel_y_centers  = tray.channel_y_centers_o
	tray_channels           = tray.channels_o
	first_hole              = tray.first_hole_i
	holes_count             = tray.holes_count_i
	frame_holes_pitch_dy    = frame.holes_pitch_dy_l
	frame_holes_dy          = frame.holes_dy_l
	rail_dx                 = west_rail.dx
	rail_wide_hole_dz       = west_rail.wide_hole_dz_l
	rail_wide_hole_diameter = west_rail.wide_hole_diameter_l

	# Define some diameters and radii:
	#tool_diameter = L(inch="1/4")
	tool_diameter = L(inch="1/8")
	tool_radius = tool_diameter/2
	post_diameter = rail_wide_hole_diameter - L(inch=0.002)
	post_radius = post_diameter/2

	# Define some X coordinates:
	zero = L()
	extra_dx = L(inch="1/4")
	x20 =  tray_dx/2
	x18 =  east_rail.e.x
	x16 =  east_rail.c.x
	x14 =  east_rail.w.x
	x10 = zero
	x6  =  west_rail.e.x
	x4  =  west_rail.c.x
	x2  =  west_rail.w.x
	x0  = -tray_dx/2

	# Define some Y coordinates:
	extra_dy = L(inch="1/4")
	y20 = tray_north_y + tool_radius
	y19 = tray_north_y
	y18 = tray_north_y - frame_holes_pitch_dy/2 + post_radius
	y16 = tray_north_y - frame_holes_pitch_dy/2
	y14 = tray_north_y - frame_holes_pitch_dy/2 - post_radius
	y10 = (tray_north_y + tray_south_y)/2
	y8  = tray_south_y + frame_holes_pitch_dy/2 + post_radius
	y6  = tray_south_y + frame_holes_pitch_dy/2
	y4  = tray_south_y + frame_holes_pitch_dy/2 - post_radius
	y1  = tray_south_y
	y0  = tray_south_y - tool_radius

	# Define some Z values:
	stock_dz = L(inch="1/2")
	extra_dz = L(inch=0.050)
	dz = stock_dz - extra_dz
	extra_top_dz = extra_dz/2
	extra_bottom_dz = zero

	# Define some Z coordinates:
	z20 = tray_base_top_z + extra_top_dz
	z18 = tray_base_top_z
	z15 = east_rail.t.z
	z10 = tray_base_top_z - dz/2
	z8  = east_rail.t.z - rail_wide_hole_dz
	z5  = east_rail.b.z
	z2  = tray_base_top_z - dz
	z0  = tray_base_top_z - dz - extra_bottom_dz

	# Create the initial block of *material*:
	material = Material("Plastic", "HDPE")
	corner1 = P(x0,  y1,  z2)
	corner2 = P(x20, y19, z18)
	tray_base.block("Block", material, color, corner1, corner2, "")

        # Mount *tray_base* into the vice, drill the tooling plate holes:
	tray_base.vice_mount("Top_Vice", "t", "n", "l",
	  extra_dx, extra_dy, extra_top_dz=extra_top_dz, extra_bottom_dz=extra_bottom_dz)
	tray_base.tooling_plate_drill("Top_Tooling_Holes", (0, 4, 8, 12, 16), (0, 4), [])

	# Remount *tray_base onto the tooling plate, mill the top surface flat, and mill
	# the exterior contour:
	tray_base.tooling_plate_mount("Top_Plate")
	tray_base.top_face("Top_Face")
	tray_base.rectangular_contour("Top_Exterior_Contour", L(inch="1/16"))

	# Now mill out the each *channel*:
	for index, channel in enumerate(tray_channels):
	    # Grab some values out of *channel* and *channel_y_centers*:
	    tape_length      = channel.tape_length
	    tape_width       = channel.tape_width
	    tape_edge_depth  = channel.tape_edge_depth
	    pocket_width     = channel.pocket_width
	    pocket_depth     = channel.pocket_width
	    pocket_offset    = channel.pocket_offset
	    channel_y_center = tray_channel_y_centers[index]
	    #print("tape_edge_depth={0:m}".format(tape_edge_depth))

	    # Compute various Y locations and Z depths:
	    tape_top_edge_y      = channel_y_center + tape_width/2
	    tape_bottom_edge_y   = channel_y_center - tape_width/2
	    pocket_top_edge_y    = tape_top_edge_y - pocket_offset + pocket_width/2
	    pocket_bottom_edge_y = tape_top_edge_y - pocket_offset - pocket_width/2

	    # Mill out a pocket that to hold the tape:
	    corner1 = P(-tape_length/2, channel_y_center - tape_width/2, z18 - tape_edge_depth)
	    corner2 = P( tape_length/2, channel_y_center + tape_width/2, z18)
	    comment = "Tray '{0}' channel {1}".format(name, index)
	    #print("corner1={0:m} corner2={1:m}".format(corner1, corner2))
	    tray_base.simple_pocket(comment, corner1, corner2, zero, "")

	    # If appropriate, mill out a pocket to hold components.
	    if pocket_width > zero:
		pocket_depth = channel.pocket_depth
		corner1 = P(-tape_length/2, pocket_bottom_edge_y, z18 - pocket_depth)
		corner2 = P( tape_length/2, pocket_top_edge_y,    z18)
		comment = "Tray '{0}' Tape channel {1}".format(name, index)
		tray_base.simple_pocket(comment, corner1, corner2, zero, "")

	# Mount the *tray_base* top side up on tooling plate:
	tray_base.vice_mount("Bottom_Vice", "b", "s", "l", extra_top_dz=extra_bottom_dz)
	tray_base.tooling_plate_drill("Bottom_Tooling_Holes", (0, 4, 8, 12, 16), (0, 4), [])
	tray_base.tooling_plate_mount("Bottom_Plate")

	# Do some common operations for each *rail*:
	fastener_diameter = "#2-56:close"
	for rail_index, rail in enumerate( [west_rail, east_rail] ):
	    x_rail_east = rail.e.x
	    x_rail_center = rail.c.x
	    x_rail_west = rail.w.x

	    # Mill out the alignment posts:
	    for y_index, y in enumerate( [y6, y16] ):
		# Drill the hole for attaching the *tray_base* to the rails:
		comment = "Hole_{0}_{1}".format(rail_index, y_index)
		start =  P(x_rail_center, y, z2)
		stop   = P(x_rail_center, y, z18)
		tray_base.hole(comment, fastener_diameter, start, stop, "t")

		# Ensure that the allignment post milled in the next operation is the right height:
		comment = "Post_Cap_{0}_{1}".format(rail_index, y_index)
		start = P(x_rail_center, y, z2)
		stop  = P(x_rail_center, y, z8)
		post_top_diameter = 1.25 * post_diameter
		tray_base.round_pocket(comment, post_top_diameter, start, stop, "")

		# Now mill out the material around the alignment post:
		comment = "Annular_Pocket_{0}_{1}".format(rail_index, y_index)
		start = P(x_rail_center, y, z2)
		stop  = P(x_rail_center, y, z15)
		outer_diameter = post_diameter + 2.40 * tool_diameter
		tray_base.round_pocket(comment,
		  outer_diameter, start, stop, "", inner_diameter=post_diameter)

	    # Mill three rail pockets between the two alignment posts:
	    for y_index, y_low_high in enumerate( [(y0, y4), (y8, y14), (y18, y20)] ):
		y_low, y_high = y_low_high
		corner1 = P(x_rail_west, y_low,  z2)
		corner2 = P(x_rail_east, y_high, z15)
		comment = "Rail_Pocket_{0}_{1}".format(rail_index, y_index)
		tray_base.simple_pocket(comment, corner1, corner2, zero, "")

class Trays(Part):
    """ *Trays*: Represents all of the trays.
    """

    def __init__(self, up, name):
	""" *Trays: Initialize the *Trays* assembly object. """

	# Standard initialization sequence:
	trays = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(trays, up, name)

	# Define each tray:
	length1 = 250.0
	length2 = 210.0
	edge_depth = 0.6
	channel8  = Channel(8.0,  length1, edge_depth,  4.55,  2.40, 1.75 +  3.50)
	channel12 = Channel(12.0, length2, edge_depth,  8.20,  6.40, 1.75 +  5.50)
	channel16 = Channel(16.0, length2, edge_depth, 12.10,  7.90, 1.75 +  7.50)
	channel24 = Channel(24.0, length2, edge_depth, 20.10, 11.90, 1.75 + 11.50)
	
	colors1 = ("lime", "green", "dark_green")
	trays.tray1_ = Tray( trays, "Tray1", colors1, 1, 2,
	  (channel8, channel8, channel8, channel8) )
	colors2 = ("pink", "red", "dark_read")
	trays.tray2_ = Tray( trays, "Tray2", colors2, 4, 2, (channel12, channel8, channel12) )
	colors3 = ("light_blue", "blue", "dark_blue")
	trays.tray3_ = Tray( trays, "Tray3", colors3, 7, 2, (channel12, channel8, channel16) )

    def construct(self):
	""" *Trays*: Construct the *Trays* asembly object."""

	pass

def main():
    ezcad = EZCAD3(0)
    small_smt = SmallSMT(None, "SmallSMT")
    small_smt.process(ezcad)

if __name__ == "__main__":
    main()
