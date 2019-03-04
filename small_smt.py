#!/usr/bin/env python

#<----------------------------------------- 100 characters --------------------------------------->|
#
# Classes are listed alphabetically.
# Methods are listed alphabetically within a class.
from EZCAD3 import *

# 10 inches of tape measure: 1070-68 = 1002 # i.e. each pix is .01in

# Rail Pitch: Notes
# Rail 1 edge values: 11 and 45
# Rail 1 center: (11+45)/2 = 28
# Rail 2 edge values: 1076 and 1110
# Rail 2 center: (1076 + 1110)/2 = 1093
# Rail center to rail center: 1093 - 28 = 1065

# Hole centers along rail:
# 155 - 37 => 1.55in - 0.37in => 29.972mm => 30mm

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
	assert isinstance(tape_length, L)
	assert isinstance(tape_edge_depth, L)
	assert isinstance(pocket_width, float)
	assert isinstance(pocket_depth, float)
	assert isinstance(pocket_offset, float)

	# Stuff values into *channel* object (i.e. *self*):
	channel = self
	channel.tape_width_l      = L(mm=tape_width)
	channel.tape_length_l     = tape_length
	channel.tape_edge_depth_l = tape_edge_depth
	channel.pocket_width_l    = L(mm=pocket_width)
	channel.pocket_depth_l    = L(mm=pocket_depth)
	channel.pocket_offset_l   = L(mm=pocket_offset)

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
	frame.rail_dx_l        = rail_dx        = L(mm=9.00)
	frame.rail_dx_pitch_l  = rail_dx_pitch  = L(mm=216.0) + rail_dx
	frame.rail_dy_l        = rail_dy        = L(mm=325.0)
	frame.holes_pitch_dy_l = holes_pitch_dy = L(mm=30.00)
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
	rail.dz_l             = dz             = L(mm=6.50)
	rail.holes_pitch_dy_l = holes_pitch_dy = L(mm=30.00)

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
	    diameter1 = L(inch=0.234)  # Letter A drill
	    #diameter2 = L(inch=0.0890)# #2-56:close
	    diameter2 = L(inch=0.1405) # #28
	    rail.hole(comment1, diameter1, start, stop1, "f")
	    rail.hole(comment2, diameter2, start, stop2, "t")
	    y += holes_pitch_dy

	# Stuff some values into *rail*:
	rail.wide_hole_dz_l         = z10 - z5
	rail.wide_hole_diameter_l   = diameter1
	rail.narrow_hole_diameter_l = diameter2

    def y_get(self, hole_index):
        """ *Rail*: Return the Y coordinate associated with *hole_index*.
	"""

	# Verify argument types:
	assert isinstance(hole_index, int)

	# Grab some *Part*'s:
	rail = self
	frame = rail.up
	
	# Compute the resulting *y* and return it:
	holes_pitch_dy = frame.holes_pitch_dy_l
	holes_dy       = frame.holes_dy_l
	y = -holes_dy/2 + hole_index * holes_pitch_dy
	return y

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

    def __init__(self, up, name, colors, hole_indices, channels):
	""" *Tray*: Initialize *Tray* assembly object.
	"""
	
	# Standard initialization sequence:
	tray = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(tray, up, name)

	# Verify additional arguments:
	assert isinstance(colors, tuple) and len(colors) == 3
	assert isinstance(hole_indices, list) or isinstance(hole_indices, tuple) 
	for hole_index in hole_indices:
	    assert isinstance(hole_index, int)
	assert isinstance(channels, tuple) and len(channels) >= 1

	# Define the tray components:
	tray.base_ = Tray_Base(tray, name + "_Base", Color(colors[0]))
	tray.lid_  = Tray_Lid(tray,  name + "_Lid",  Color(colors[1]))
	tray.cap_  = Tray_Cap(tray,  name + "_Cap",  Color(colors[2]))

	# Save the additional arguments into *tray*:
	tray.hole_indices_o = hole_indices
	tray.channels_o     = channels

	# Initialize the various center and width lists:
	channels_size = len(channels)
	zero                     = L()
	tray.channel_y_centers_o = [zero] * channels_size
	tray.lid_y_centers_o     = [zero] * channels_size
	tray.lid_y_widths_o      = [zero] * channels_size
	tray.gap_y_centers_o     = [zero] * (channels_size + 1)

    def construct(self):
	""" *Tray*: Constructs the *Tray* object. """

	# Grab some values from *tray* (i.e. *self*):
	tray              = self
	channels          = tray.channels_o
	channel_y_centers = tray.channel_y_centers_o
	hole_indices      = tray.hole_indices_o
	gap_y_centers     = tray.gap_y_centers_o
	lid_y_centers     = tray.lid_y_centers_o
	lid_y_widths      = tray.lid_y_widths_o

	# Grab some useful *Part*'s:
	trays     = tray.up
	small_smt = trays.up
	frame     = small_smt.frame_
	rail      = frame.west_rail_ # Either the east or west rail would work here.
	
	# Grab some useful values from *frame*:
	rail_dx        = frame.rail_dx_l
	rail_dx_pitch  = frame.rail_dx_pitch_l
	holes_pitch_dy = frame.holes_pitch_dy_l
	holes_dy       = frame.holes_dy_l

	# Determine the *first_hole_index* and *last_hole_index* and save into *tray*:
	tray.high_hole_index_i = high_hole_index = max(hole_indices)
	tray.low_hole_index_i  = low_hole_index  = min(hole_indices)
	tray.holes_span_i      = holes_span      = high_hole_index - low_hole_index

	# Define some Y and Z values:
	zero = L()
	tray.shave_dy_l   = shave_dy   = L(inch=0.005)
	tray.dx_l         = dx         = rail_dx_pitch + 3 * rail_dx 
	tray.north_y_l    = north_y    = rail.y_get(high_hole_index) + holes_pitch_dy/2 - shave_dy
	tray.south_y_l    = south_y    = rail.y_get(low_hole_index)  - holes_pitch_dy/2 + shave_dy
	tray.dy_l         = dy         = north_y - south_y
	tray.base_top_z_l = base_top_z = zero   # Top surface of *Tray_Base*
	#print("north_y={0:m} south_y={1:m}".format(north_y, south_y))

	# Compute *total_channels_width*:
	total_channels_width = zero
	for channel_index, channel in enumerate(channels):
	    total_channels_width += channel.tape_width_l
	    #print("[{0}]: total_channels_width={1:m}".format(channel_index, total_channels_width))
	#print("total_channels_width={0:m}".format(total_channels_width))

	# *edge_gap_dy* is the gap on the edge of part and is fixed.  *center_gap_dy* is
        # is the gap between the channels and is computed below:
	zero = L()
	channels_size = len(channels)
	if channels_size == 1:
	    edge_gap_dy = (dy - total_channels_width) / 2
	    center_gap_dy = None
	else:
	    edge_gap_dy = L(mm=5.00)
	    center_left_over = dy - 2 * edge_gap_dy - total_channels_width
	    assert center_left_over > zero
	    center_gap_dy = center_left_over / (channels_size - 1)
	    #print("dy={0:m} center_left_over={1:m} center_gap_dy={2:m}".
	    #  format(dy, center_left_over, center_gap_dy))

	# Now compute in the center line for each channel and stuff into *channel_y_centers*.
	# Also compute the center line for each gap in *gap_y_centers*:
	y = south_y
	channel_y_centers = tray.channel_y_centers_o
	for channel_index, channel in enumerate(channels):
	    # The first time through *gap_dy* is *edge_gap_dy* and all other times *center_gap_dy*:
	    gap_dy = edge_gap_dy if channel_index == 0 else center_gap_dy

	    # Compute the next value for *gap_y_centers*:
	    gap_y_centers[channel_index] = y + gap_dy/2
	    y += gap_dy
	
	    # Compute the next value for *channel_y_centers*:
	    tape_width = channel.tape_width_l
	    channel_y_centers[channel_index] = y + tape_width/2
	    y += tape_width

	    #print("[{0}] y={1:m} ".format(channel_index, y))

	# Fill in the last value for *gap_y_centers*:
	gap_y_centers[channels_size] = north_y - edge_gap_dy/2

	# Now figure out the lid boundaries:
	lid_y_centers = tray.lid_y_centers_o
	lid_y_widths  = tray.lid_y_widths_o
	for channel_index, channel in enumerate(channels):
	    tape_width   = channel.tape_width_l
	    channel_y_center = channel_y_centers[channel_index]
	    tape_north_y = channel_y_center + tape_width / 2
	    tape_south_y = channel_y_center - tape_width / 2

	    # Compute the default *lid_y_center* and *lid_y_width* for when there is no pocket.
	    pocket_width  = channel.pocket_width_l
	    # 1.75mm = north edge to holes center
	    # 1.50mm/2 = holes center to bottom hole edge
	    # 1.50mm/4 = bottom hole edge to estimated top cover tape edge
	    lid_y_width = tape_width - L(mm=1.75 + 1.50/2 + 1.50/4)
	    lid_y_center = tape_south_y + lid_y_width / 2

	    # If we have a pocket, we compute a different *lid_y_center* and *lid_y_width*:
	    pocket_width  = channel.pocket_width_l
	    if pocket_width > zero:
		# We have a pocket, so we compute a different lid_y
		pocket_offset = channel.pocket_offset_l
		lid_center_y = tape_north_y - pocket_offset
		lid_y_width = pocket_width

	    # Stuff values into *lid_y_centers* and *lid_y_widths*:
	    lid_y_centers[channel_index] = lid_y_center
	    lid_y_widths[channel_index]  = lid_y_width
	    #print("[{0}]:lid_y_center={1:m} lid_y_width={2:m}".
	    #  format(channel_index, lid_y_center, lid_y_width))

    def base_mount_holes_drill(self, part, diameter):
	""" *Tray*: Drill the mount holes for the *Channel* object.

	The arguments are:
	* *self* (i.e. *Tray*) : The *Tray* object to fetch the *Channel* objects from.
	* *part* (*Part*): The part to drill the holes into.
	* *diameter* (*L*): The diameter of the holes.
	"""

	# Verify argument types:
	assert isinstance(part, Tray_Base) \
	  or isinstance(part, Tray_Lid) or isinstance(part, Tray_Cap)
	assert isinstance(diameter, L) or isinstance(diameter, str)
	
	# Grap some *Part*'s from *tray*:
	tray      = self
	trays     = tray.up
	small_smt = trays.up
	frame     = small_smt.frame_
	east_rail = frame.east_rail_
	west_rail = frame.west_rail_

	# Grab some values from *tray*:
	tray_hole_indices    = tray.hole_indices_o
	tray_high_hole_index = tray.high_hole_index_i
	tray_low_hole_index  = tray.low_hole_index_i
	part_name = part.name_get()

	# Get the top and bottom *part*:
	z10 = part.t.z
	z0  = part.b.z

	# Do some common operations for each *rail*:
	for rail_index, rail in enumerate( [west_rail, east_rail] ):
	    x = rail.c.x
	    for index, hole_index in enumerate(tray_hole_indices):
		y = rail.y_get(hole_index)
		start = P(x, y, z10)
		stop  = P(x, y, z0)
		comment = "{0}_{1}_{2}".format(part_name, rail_index, index)
		part.hole(comment, diameter, start, stop, "t")

    def cap_mount_holes_drill(self, part, diameter):
	""" *Tray*: Drill the holes to attach the *Tray_Cap* to the *Tray_Base*.

	The arguments are:
	* *self* (i.e. *Tray*) : The *Tray* object to fetch the *Channel* objects from.
	* *_part* (*Part*): The part to drill the holes into.
	* *diameter* (*L*): The diameter of the holes.
	"""

	# Verify argument types:
	assert isinstance(part, Tray_Base) \
	  or isinstance(part, Tray_Lid) or isinstance(part, Tray_Cap)
	assert isinstance(diameter, L) or isinstance(diameter, str)
	
	# Grab some *Part*'s from *tray* (i.e. *self*):
	tray                 = self
	trays                = tray.up
	small_smt            = trays.up
	frame                = small_smt.frame_
	rail                 = frame.west_rail_ # Either *Rail* will do

	# Grap some values from *tray* and *part*:
	tray_name            = tray.name_get()
	tray_high_hole_index = tray.high_hole_index_i
	tray_low_hole_index  = tray.low_hole_index_i
	part_name            = part.name_get()

	# Define some X/Y/Z coordinates:
	x_offset = L(mm=5.00)
	x10 = part.e.x - x_offset
	x0  = part.w.x + x_offset

	if tray_low_hole_index == tray_high_hole_index:
	    y = rail.y_get(tray_low_hole_index)
	    y_offset = L(mm=10.00)
	    y10 = y + y_offset
	    y0  = y - y_offset
	else:
	    y10 = rail.y_get(tray_low_hole_index)
	    y0  = rail.y_get(tray_high_hole_index)

	z10 = part.t.z
	z0  = part.b.z

	# Drill holes for each *channel*:
	for x_index, x in enumerate( [x0, x10] ):
	    for y_index, y in enumerate( [y0, y10] ):
		start = P(x, y, z10)
		stop  = P(x, y, z0)
		comment = "Hole_{0}_{1}_{2}".format(part_name, x_index, y_index)
		part.hole(comment, diameter, start, stop, "t")

    def lid_cap_holes_drill(self, lid_cap_part, diameter):
	""" *Tray*: Drill the lid/cap holes that attach the "filler" pieces to the cap:

	The arguments are:
	* *self* (i.e. *Tray*) : The *Tray* object to fetch the *Channel* objects from.
	* *lid_cap_part* (*Part*): The part to drill the holes into.
	* *diameter* (*L*): The diameter of the holes.
	"""

	# Verify argument types:
	assert isinstance(lid_cap_part, Tray_Lid) or isinstance(lid_cap_part, Tray_Cap)
	assert isinstance(diameter, L) or isinstance(diameter, str)
	
	# Grap some values from *tray* and *lid_cap_part*:
	tray                = self
	trays               = tray.up
	channels            = tray.channels_o
	channel_y_centers   = tray.channel_y_centers_o
	lid_y_centers       = tray.lid_y_centers_o
	lid_y_widths        = tray.lid_y_widths_o
	tray_name           = tray.name_get()

	trays_normal_length = trays.normal_length_l
	lid_cap_part_name   = lid_cap_part.name_get()

	# Figure the top and bottom of *part*:
	z10 = lid_cap_part.t.z
	z0  = lid_cap_part.b.z

	# Drill holes for each *channel*:
	for channel_index, channel in enumerate(channels):
	    # Define some X coordinate values:
	    zero = L()
	    x_offset = L(mm=10.00)
	    x10 =  trays_normal_length/2 - x_offset
	    x5  =  zero
	    x0  = -trays_normal_length/2 + x_offset

	    # Grab some values from *channel*:
	    lid_y_center = lid_y_centers[channel_index]	    
	    lid_y_width  = lid_y_widths[channel_index]	    
	    y10 = lid_y_center + lid_y_width/2
	    y5  = lid_y_center
	    y0  = lid_y_center + lid_y_width/2

	    # Drill 3 holes across each *channel*:
	    for x_index, x in enumerate( [x0, x5, x10] ):
		start = P(x, y5, z10)
		stop  = P(x, y5, z0)
                comment = "Hole_{0}_{1}_{2}".format(lid_cap_part_name, channel_index, x_index)
		#tracing = 3 if (isinstance(lid_cap_part, Tray_Lid) and channel_index == 0 and
		#  x_index == 0 and tray_name == "Tray1" ) else -1000000
		lid_cap_part.hole(comment, diameter, start, stop, "t") #, tracing=tracing)

    def lid_mount_holes_drill(self, part, diameter):
	""" *Tray*: Drill the holes to mount the *Cap_Lid* to the *Cap_Base*.

	The arguments are:
	* *self* (i.e. *Tray*) : The *Tray* object to fetch the *Channel* objects from.
	* *part* (*Part*): The part to drill the holes into.
	* *diameter* (*L*): The diameter of the holes.
	"""

	# Verify argument types:
	assert isinstance(part, Tray_Base) \
	  or isinstance(part, Tray_Lid) or isinstance(part, Tray_Cap)
	assert isinstance(diameter, L) or isinstance(diameter, str)
	
	#print("=>Tray.lid_mount_holes('{0}', '{1}', {2})".format(
	#  self.name_get(), part.name_get(), diameter))

	# Grap some values from *tray* and *lid_cap_part*:
	tray                = self
	trays               = tray.up
	channels            = tray.channels_o
	gap_y_centers       = tray.gap_y_centers_o
	tray_name           = tray.name_get()
	trays_normal_length = trays.normal_length_l
	part_name           = part.name_get()

	# Figure the top and bottom of *part*:
	z10 = part.t.z
	z0  = part.b.z

	# Drill holes for each *channel*:
	for gap_index, gap_y_center in enumerate(gap_y_centers):
	    # Define some X coordinate values:
	    zero = L()
	    x_offset = L(mm=10.00)
	    x10 =  trays_normal_length/2 - x_offset
	    x5  =  zero
	    x0  = -trays_normal_length/2 + x_offset

	    # Grab some values from *channel*:
	    y = gap_y_center

	    # Drill 3 holes across each *channel*:
	    for x_index, x in enumerate( [x0, x5, x10] ):
		start = P(x, y, z10)
		stop  = P(x, y, z0)
                comment = "Hole_{0}_{1}_{2}".format(part_name, gap_index, x_index)
		#tracing = 3 if (part.name_get() == "Tray889_Base" and
		#  isinstance(diameter, str) and diameter == "#2-56:thread" and
		#  gap_index == 0 and x_index< == 0) else -1000000
		part.hole(comment, diameter, start, stop, "t") #, tracing=tracing)

    def pin_holes_drill(self, part, diameter):
	""" *Tray*: Drill the pin alignment holes for the *Channel* object.

	The arguments are:
	* *self* (i.e. *Tray*) : The *Tray* object to fetch the *Channel* objects from.
	* *part* (*Part*): The part to drill the holes into.
	* *diameter* (*L*): The diameter of the holes.
	"""

	# Verify argument types:
	assert isinstance(part, Tray_Base) or isinstance(part, Tray_Lid)
	assert isinstance(diameter, L) or isinstance(diameter, str)
	
	# Grap some values from *tray* and *part*:
	tray                = self
	trays               = tray.up
	channels            = tray.channels_o
	channel_y_centers   = tray.channel_y_centers_o
	tray_name           = tray.name_get()
	trays_normal_length = trays.normal_length_l
	part_name           = part.name_get()

	# The pin to pin spacing is 4mm.  We want to have pin alignment holes on a regular
        # basis.  5 * 4 = 20, so every 5 pin hole has an alignment hole:
	pin_pitch_dx = L(mm=20.00)
	pins_span = int(trays_normal_length / pin_pitch_dx)
	pins_dx = pins_span * pin_pitch_dx

	# Figure the top and bottom of *part*:
	part_top_z    = part.t.z
	part_bottom_z = part.b.z

	# Drill holes for each *channel*:
	for channel_index, channel in enumerate(channels):
	    # Grab some values from *channel*:
	    channel_y_center = channel_y_centers[channel_index]
	    tape_edge_depth = channel.tape_edge_depth_l
	    tape_width = channel.tape_width_l
	    tape_north_y = channel_y_center + tape_width/2
	    y  = tape_north_y - L(mm=1.75)

	    # The bottom of the drill depends upon whether *part* is a *Tray_Base*
            # or not.  For a *Tray_Base*, we only want to go 2.5mm below the *tape_edge_depth*.
            # Otherwise, we want to go all the way through the lid:
	    if isinstance(part, Tray_Base):
		tape_bottom_depth_z = part_top_z - tape_edge_depth
		z_drill_bottom = part_top_z - tape_edge_depth - L(mm=2.50)
		hole_flags = "p"
	    elif isinstance(part, Tray_Lid):
		z_drill_bottom = part_bottom_z
		hole_flags = "t"
	    else:
		assert False, "This should not be possible"

	    # Drill 3 holes across each *channel*:
	    for x_index, x in enumerate(range(pins_span + 1)):
		x = -pins_dx/2 + x_index * pin_pitch_dx
		start = P(x, y, part_top_z)
		stop  = P(x, y, z_drill_bottom)
                comment = "Hole_{0}_{1}_{2}".format(part_name, channel_index, x_index)
		#tracing = 3 if (isinstance(lid_cap_part, Tray_Lid) and channel_index == 0 and
		#  x_index == 0 and tray_name == "Tray1" ) else -1000000

		# Take the drill down to *stop* in tip mode with `"p"` flag:
		part.hole(comment, diameter, start, stop, hole_flags) #, tracing=tracing)

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
	tray_high_hole_index    = tray.high_hole_index_i
	tray_hole_indices       = tray.hole_indices_o
	tray_holes_span         = tray.holes_span_i
	tray_low_hole_index     = tray.low_hole_index_i
	frame_holes_pitch_dy    = frame.holes_pitch_dy_l
	frame_holes_dy          = frame.holes_dy_l
	rail_dx                 = west_rail.dx
	rail_wide_hole_dz       = west_rail.wide_hole_dz_l
	rail_wide_hole_diameter = west_rail.wide_hole_diameter_l

	# Define some diameters and radii:
	#tool_diameter = L(inch="1/4")
	tool_diameter = L(inch="3/16")
	#tool_diameter = L(inch="1/8")
	tool_radius = tool_diameter/2
	post_diameter = rail_wide_hole_diameter - L(inch=0.001)
	post_radius = post_diameter/2

	# Define some X coordinates:
	zero = L()
	extra_dx = L(inch="1/4")
	x20 =  tray_dx/2
	x18 =  east_rail.e.x
	x16 =  east_rail.c.x
	x14 =  east_rail.w.x
	x10 =  zero
	x6  =  west_rail.e.x
	x4  =  west_rail.c.x
	x2  =  west_rail.w.x
	x0  = -tray_dx/2

	# Define some Y coordinates:
	north_hole_y = west_rail.y_get(tray_high_hole_index)
	south_hole_y = west_rail.y_get(tray_low_hole_index)
	extra_dy = L(inch="1/4")
	y20 = tray_north_y + tool_radius
	y19 = tray_north_y
	y10 = (tray_north_y + tray_south_y)/2
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

	tooling_plate_columns = (0, 4, 7, 9, 12, 16)
	if tray_holes_span == 0:
	    tooling_plate_rows = (0, 2)
	elif tray_holes_span == 1:
	    tooling_plate_rows = (0, 4)
	elif tray_holes_span == 2:
            tooling_plate_rows = (0, 6)

        # Mount *tray_base* into the vice, drill the tooling plate holes:
	tray_base.vice_mount("Top_Vice", "t", "n", "l",
	  extra_dx, extra_dy, extra_top_dz=extra_top_dz, extra_bottom_dz=extra_bottom_dz)
	tray_base.tooling_plate_drill("Top_Tooling_Holes",
	  tooling_plate_columns, tooling_plate_rows, [])

	# Remount *tray_base onto the tooling plate, mill the top surface flat, and mill
	# the exterior contour:
	tray_base.tooling_plate_mount("Top_Plate")
	tray_base.top_face("Top_Face")
	tray_base.rectangular_contour("Top_Exterior_Contour", L(inch="1/16"))

	# Now mill out the each *channel*:
	for index, channel in enumerate(tray_channels):
	    # Grab some values out of *channel* and *channel_y_centers*:
	    tape_width       = channel.tape_width_l
	    tape_length      = channel.tape_length_l
	    tape_edge_depth  = channel.tape_edge_depth_l
	    pocket_width     = channel.pocket_width_l
	    pocket_depth     = channel.pocket_width_l
	    pocket_offset    = channel.pocket_offset_l
	    channel_y_center = tray_channel_y_centers[index]
	    pocket_depth     = channel.pocket_depth_l
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
		corner1 = P(-tape_length/2, pocket_bottom_edge_y, z18 - pocket_depth)
		corner2 = P( tape_length/2, pocket_top_edge_y,    z18)
		comment = "Tray '{0}' Tape channel {1}".format(name, index)
		tray_base.simple_pocket(comment, corner1, corner2, zero, "")

	# Make sure that the milling is done before doing the drilling:
	tray_base.cnc_fence()

	# Drill the pin alignment holes.  Note that "#0-80:close" is a #52 drill which is
        # 0.635 inches in diameter  which is just a little bit bigger than 1.5mm=(.0590in)
        # which is the pin hole diameter.  The #52 drill is premounted in a tool holder
        # for Wayne's mill:
	tray.pin_holes_drill(tray_base, "#0-80:close")

	# Drill the cap and lid attach holes:
	tray.cap_mount_holes_drill(tray_base, "#2-56:thread")
	tray.lid_mount_holes_drill(tray_base, "#2-56:thread")
	# Note: the base mount holes are drilled after we flip the *tray_base* over:

	# Now mount the *tray_base* bottom side up on tooling plate:
	tray_base.vice_mount("Bottom_Vice", "b", "s", "l", extra_top_dz=extra_bottom_dz)
	tray_base.tooling_plate_drill("Bottom_Tooling_Holes",
	  tooling_plate_columns, tooling_plate_rows, [])
	tray_base.tooling_plate_mount("Bottom_Plate")

	# Do some common operations for each *rail*:
	fastener_diameter = "#2-56:close"
	for rail_index, rail in enumerate( [west_rail, east_rail] ):
	    x_rail_east = rail.e.x
	    x_rail_center = rail.c.x
	    x_rail_west = rail.w.x

	    # Mill out the alignment posts:
	    for index, hole_index in enumerate(tray_hole_indices):
		y = west_rail.y_get(hole_index)

		# Drill the hole for attaching the *tray_base* to the rails:
		comment = "Hole_{0}_{1}".format(rail_index, hole_index)
		start =  P(x_rail_center, y, z2)
		stop   = P(x_rail_center, y, z18)
		tray_base.hole(comment, fastener_diameter, start, stop, "t")

		# Ensure that the allignment post milled in the next operation is the right height:
		comment = "Post_Cap_{0}_{1}".format(rail_index, hole_index)
		start = P(x_rail_center, y, z2)
		stop  = P(x_rail_center, y, z8)
		post_top_diameter = 1.25 * post_diameter
		tray_base.round_pocket(comment, post_top_diameter, start, stop, "")

		# Now mill out the material around the alignment post:
		comment = "Annular_Pocket_{0}_{1}".format(rail_index, hole_index)
		start = P(x_rail_center, y, z2)
		stop  = P(x_rail_center, y, z15)
		# 2.06? Yeah it's a kludge.  This deals with the fact that the *round_pocket*
		# operation needs some extra space for the "spring pass" which is currently .005in.
		# If we don't multiply by 2.06, we get the laser instead of the 3/16 end-mill:
		outer_diameter = post_diameter + 2.06 * tool_diameter
		#tracing = 3 if rail_index == 0 and index == 0 else -1000000
		tray_base.round_pocket(comment,
		  outer_diameter, start, stop, "", inner_diameter=post_diameter) #, tracing=tracing)

	    # Mill rail pockets between the alignment posts and at each end:
	    hole_indices = tray_hole_indices
	    hole_indices_size = len(hole_indices)
	    rail_adjust = L(inch=0.002)
	    for index in range(hole_indices_size + 1):
		# If *index* is 0, we do the bottom end rail pocket.
                # If *index* is *hole_indices_size*, we do the top end rail pocket.
		# In all other cases we do the pockts between two posts:
		y_low  = ( tray_south_y - tool_radius if index == 0
		  else rail.y_get(hole_indices[index - 1]) + post_radius)
		y_high = ( tray_north_y + tool_radius if index == hole_indices_size
		  else rail.y_get(hole_indices[index])     - post_radius)
		#print("[{0}, {1}]: y_low={2:m} y_high={3:m}".
		#  format(rail_index, index, y_low, y_high))
		corner1 = P(x_rail_west - rail_adjust, y_low,  z2)
		corner2 = P(x_rail_east + rail_adjust, y_high, z15)
		comment = "Rail_Pocket_{0}_{1}".format(rail_index, index)
		tray_base.simple_pocket(comment, corner1, corner2, tool_radius, "")

class Tray_Lid(Part):
    """ *Tray_Lid*: Represent the lid that holes the cut tape in place in the *Tray_Base*.
    """

    def __init__(self, up, name, color):
	""" *Tray_Lid: Initialize the *Tray_Lid* object (i.e. *self*.) """

	# Standard initialization sequence:
	tray_lid = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(color, Color)
	Part.__init__(tray_lid, up, name)
	tray_lid.color_o = color

    def construct(self):
        """ *Tray_Lid*: construct the *Tray_Lid* object (i.e. *self*.)
	"""

	# Grab some *Part*'s from the *Tray_Lid* object (i.e. *self*):
	tray_lid = self
	tray     = tray_lid.up
	trays    = tray.up

	# Grap some values from the *Part*'s:
	trays_normal_length    = trays.normal_length_l
	trays_debug            = trays.debug_b
	tray_channels          = tray.channels_o
	tray_channel_y_centers = tray.channel_y_centers_o
	tray_dx                = tray.dx_l	
	tray_dy                = tray.dy_l
	tray_north_y           = tray.north_y_l
	tray_south_y           = tray.south_y_l
	tray_base_top_z        = tray.base_top_z_l
	color                  = tray_lid.color_o

	tray_lid.dz_l = dz = L(inch="1/8")

	x10 =  tray_dx/2
	x8  =  tray_dx/4
	x2  = -tray_dx/4
	x0  = -tray_dx/2

	y_offset = L(inch="1/4")
	y10 = tray_north_y
	y8  = tray_north_y - y_offset
	y2  = tray_south_y + y_offset
	y0  = tray_south_y

	z10 = tray_base_top_z + dz
	z0  = tray_base_top_z

	material = Material("Plastic", "HDPE")
	name = tray_lid.name_get()
	comment = "Tray_Lid_{0}".format(name)
	corner1 = P(x0,  y0,  z0)
	corner2 = P(x10, y10, z10)
	zero = L()
	tray_lid.block(comment, material, color, corner1, corner2, "")
	tray_lid.vice_mount("Top_Vice", "t", "n", "", zero, zero)
	tray_lid.tool_prefer("Laser_007")
	tray_lid.rectangular_contour("Rectangular Contour", L(inch="1/16"))

	# Now mill out the each *channel*:
	for index, channel in enumerate(tray_channels):
	    # Grab some values out of *channel* and *channel_y_centers*:
	    tape_width       = channel.tape_width_l
	    tape_edge_depth  = channel.tape_edge_depth_l
	    pocket_width     = channel.pocket_width_l
	    pocket_depth     = channel.pocket_width_l
	    pocket_offset    = channel.pocket_offset_l
	    channel_y_center = tray_channel_y_centers[index]
	    #print("tape_edge_depth={0:m}".format(tape_edge_depth))

	    # Compute various Y locations and Z depths:
	    tape_top_edge_y      = channel_y_center + tape_width/2
	    tape_bottom_edge_y   = channel_y_center - tape_width/2
	    pocket_top_edge_y    = tape_top_edge_y - pocket_offset + pocket_width/2
	    pocket_bottom_edge_y = tape_top_edge_y - pocket_offset - pocket_width/2

	    # Mill out a pocket that to hold the tape:
	    corner1 = P(-trays_normal_length/2, pocket_bottom_edge_y, z0)
	    corner2 = P( trays_normal_length/2, pocket_top_edge_y,    z10)
	    comment = "Tray '{0}' channel {1}".format(name, index)
	    #print("corner1={0:m} corner2={1:m}".format(corner1, corner2))
	    tray_lid.simple_pocket(comment, corner1, corner2, zero, "")

	# "#0-80:thread" is same as the hole drill in the *Tray_Base*:
	tray.pin_holes_drill(tray_lid, "#0-80:thread")

	tray.lid_cap_holes_drill(tray_lid,    "#2-56:thread")

	# Drill the various mounting hiles holes for the *tray_lid*:
	tray.base_mount_holes_drill(tray_lid, "#2-56:pan_head")
	tray.lid_mount_holes_drill(tray_lid,  "#2-56:close")
	tray.cap_mount_holes_drill(tray_lid,  "#2-56:close")

	# Perform an requested debug visualization:
	if trays_debug:
	    corner_radius = L(inch="1/4")
	    corner1 = P(x2, y2, z0)
	    corner2 = P(x8, y8, z10)
	    tray_lid.simple_pocket("debug pocket", corner1, corner2, corner_radius, "t")

class Tray_Cap(Part):
    """ *Tray_Cap*: Represents the cap covers everything when a *Tray* is being moved around.
    """

    def __init__(self, up, name, color):
	""" *Tray_Cap: Initialize the *Tray_Cap* object (i.e. *self*.) """

	# Standard initialization sequence:
	tray_cap = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(color, Color)
	Part.__init__(tray_cap, up, name)
	tray_cap.color_o = color

    def construct(self):
        """ *Tray_Cap*: construct the *Tray_Cap* object (i.e. *self*.)
	"""

	# Grab some *Part*'s from the *Tray_Cap* object (i.e. *self*):
	tray_cap = self
	tray     = tray_cap.up
	trays    = tray.up
	tray_lid = tray.lid_

	# Grap some values from the *Part*'s:
	trays_normal_length = trays.normal_length_l
	trays_debug         = trays.debug_b
	tray_channels       = tray.channels_o
	tray_dx             = tray.dx_l	
	tray_dy             = tray.dy_l
	tray_north_y        = tray.north_y_l
	tray_south_y        = tray.south_y_l
	tray_base_top_z     = tray.base_top_z_l
	color               = tray_cap.color_o

	# Define some X/Y/Z coordinates:
	tray_lid.dz_l = dz = L(inch="1/8")
	x10 =  tray_dx/2
	x8  =  tray_dx/3
	x2  = -tray_dx/3
	x0  = -tray_dx/2
	y_offset = L(inch="1/4")
	y10 = tray_north_y
	y8  = tray_north_y - y_offset
	y2  = tray_south_y + y_offset
	y0  = tray_south_y
	z10 = tray_lid.t.z + dz
	z0  = tray_lid.t.z

	# Construct the original block out of *material*:
	material = Material("Plastic", "HDPE")
	name = tray_cap.name_get()
	comment = "Tray_Lid_{0}".format(name)
	corner1 = P(x0,  y0,  z0)
	corner2 = P(x10, y10, z10)
	zero = L()
	tray_cap.tool_prefer("Laser_007")
	tray_cap.block(comment, material, color, corner1, corner2, "")
	tray_cap.vice_mount("Top_Vice", "t", "n", "", zero, zero)
	tray_cap.rectangular_contour("Rectangular Contour", L(inch="1/16"))

	# Drillthe holes to attach the "filler" pieces to the bottom of *tray_cap*:
	tray.lid_cap_holes_drill(tray_cap,    "#2-56:close")

	# Drill the various mounting holes into *tray_cap*:
	tray.base_mount_holes_drill(tray_cap, "#2-56:pan_head")
	tray.lid_mount_holes_drill(tray_cap,  "#2-56:pan_head")
	tray.cap_mount_holes_drill(tray_cap,  "#2-56:close")

	# Perform any requested debug visualization:
	if trays_debug:
	    corner_radius = L(inch="1/4")
	    corner1 = P(x2, y2, z0)
	    corner2 = P(x8, y8, z10)
	    tray_cap.simple_pocket("Debug pocket", corner1, corner2, corner_radius, "t")

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

	# Some *trays* values:
	trays.full_length_l   = full_length   = L(mm=270.0)
	trays.normal_length_l = normal_length = L(mm=209.0)
	trays.debug_b         = debug         = True

	# Define all of the channels:
	edge_depth = L(mm=0.6)
	channel8f  = Channel(8.0, full_length,    L(mm=0.64),  4.55,  2.40, 1.75 +  3.50)
	channel8n  = Channel(8.0, normal_length,  L(mm=0.64),  4.55,  2.40, 1.75 +  3.50)
	channel12n = Channel(12.0, normal_length, edge_depth,  8.20,  6.40, 1.75 +  5.50)
	channel12f = Channel(12.0, full_length,   edge_depth,  8.20,  6.40, 1.75 +  5.50)
	channel16 = Channel(16.0, normal_length,  edge_depth, 12.10,  7.90, 1.75 +  7.50)
	channel24 = Channel(24.0, normal_length,  edge_depth, 20.10, 11.90, 1.75 + 11.50)

	# Finally define all the *Tray* *Part*'s:
	colors1 = ("lime", "green", "dark_green")
	#trays.tray1_ = Tray( trays, "Tray889", colors1, (1, 3),
	#  (channel8n, channel8f, channel8f, channel8f, channel8n) )
	colors2 = ("pink", "red", "dark_red")
	#trays.tray2_ = Tray( trays, "Tray2222", colors2, (4, 6),
	#  (channel12n, channel12f, channel12f, channel12n) )
	colors3 = ("light_blue", "blue", "dark_blue")
	#trays.tray3_ = Tray( trays, "Tray666", colors3, (7, 9),
	# (channel16, channel16, channel16) )
	colors4 = ("aqua", "medium_purple", "purple")
	trays.tray4_ = Tray( trays, "Tray211", colors4, (10,),  (channel12n, ) )

    def construct(self):
	""" *Trays*: Construct the *Trays* asembly object."""

	# This is a pure assembly, so there is nothing to construct:
	pass

def main():
    ezcad = EZCAD3(0)
    small_smt = SmallSMT(None, "SmallSMT")
    small_smt.process(ezcad)

if __name__ == "__main__":
    main()
