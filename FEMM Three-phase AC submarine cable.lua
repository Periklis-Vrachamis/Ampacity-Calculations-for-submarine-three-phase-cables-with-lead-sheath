-- =========================
-- Power Core Cable (Heat Flow FEMM) https://www.femm.info/wiki/HomePage
-- Uses Lua 4.0 math library per: https://www.lua.org/manual/4.0/manual.html#mathlib
-- Geometry in millimeters
-- =========================

-- -------------------------
-- Given cable geometry parameters (mm) 
-- -------------------------
conductor_diameter = 30              -- Conductor diameter
conductor_tape_diameter = 30.5       -- External diameter of conductor tape
conductor_screen_outer_diameter = 32.5 -- External diameter of sc conductor screen
insulation_outer_diameter = 48.5     -- External diameter of insulation
insulation_screen_outer_diameter = 50.3 -- External diameter of sc insulation screen
water_blocking_outer_diameter = 52.5 -- External diameter of sc water blocking tape
lead_sheath_outer_diameter = 57.1    -- External diameter of metallic (lead) sheath
pe_sheath_outer_diameter = 62.1      -- External diameter of sc PE sheath
filler_diameter = 134.9                  -- External diameter of fillers (cable after assembling)
cable_core_binder_diameter = 135.3     -- External diameter of cable core binder
armour_bedding_diameter = 137.3          -- External diameter of armour bedding
armour_outer_diameter = 149.3            -- External diameter of the armour
cable_outer_diameter = 155.3             -- External diameter of outer covering

-- Armour wires Specificatons
material_armour = "Galvanized Stainless Steel" -- Set as "Galvanized Stainless Steel" OR "Galvanized Steel"
number_armour_wires = 71 -- Number of Armour wires
diameter_armour_wires = 6  -- Diameter of Armour wires (mm)
armour_lay_length = 1785 -- Lay length of the Armour (mm)

-- System parameters 
Current = 974                          -- Steady state current (A)
burial_depth = 1000                    -- Burial depth to cable center (mm) 
seabed_temperature = 15.0              -- Seabed soil temperature (°C) 
seabed_thermal_resistivity = 0.7       -- Soil thermal resistivity (K·m/W)
number_cables = 1                      -- Number of cables. Maximum 3 
axial_separation  = 0                  -- Axial separation between cables (mm)
system_voltage = 30000                -- System voltage (V)
system_frequency = 50                  -- System frequency (Hz)
R_dc_20 = 0.0283e-3                    -- Ω/m 
alpha = 3.93e-3                        -- Temperature coefficient based on IEC 60287-1-1 Table 1
skin_effect_factor = 1                 -- Based on IEC 60287-1-1 Table 2
proximity_effect_factor = 1            -- Based on IEC 60287-1-1 Table 2
epsilon = 2.5                          -- Relative permittivity of XLPE based on IEC 60287-1-1 Table 3
tan_delta = 4e-3                       -- Loss factor of XLPE (filled and unfilled) based on IEC 60287-1-1 Table 3
resistivity_sheath_20 = 21.4e-8        -- Ω·m at 20°C (electrical resistivity) based on IEC 60287-1-1 Table 1
alpha_sheath = 4.0e-3                  -- Temperature coefficient based on IEC 60287-1-1 Table 1
resistivity_armour_20 = 13.8e-8        -- Ω·m at 20°C (electrical resistivity) based on IEC 60287-1-1 Table 1
alpha_armour = 4.5e-3                  -- Temperature coefficient based on IEC 60287-1-1 Table 1

-- Core layup parameters - 
core_lay_length = 2152                 -- Lay length of core stranding (mm)
phase_spacing = pe_sheath_outer_diameter    -- Distance between conductor axes 

-- Thermal resistivities (K·m/W) 
rho_conductor = 1.0/385.0              -- Copper (assumed)
rho_conductor_tape = 6                 -- Binding tape (CIGRE TB 880 Note 2 - Table 12-11)
rho_conductor_screen = 2.5             -- SC conductor screen (IEC 60287-2-1 Table 1 and CIGRE TB 880 Table 2-3)
rho_xlpe_insulation = 3.5              -- XLPE insulation (IEC 60287-2-1 Table 1)
rho_insulation_screen = 2.5            -- SC insulation screen (IEC 60287-2-1 Table 1 and CIGRE TB 880 Table 2-3)  
rho_water_blocking = 12.0              -- Water blocking tapes (CIGRE TB 880 Table 2-3)
rho_lead_sheath = 1.0/46.7             -- Lead sheath (assumed)
rho_pe_sheath = 2.5                    -- SC PE sheath (CIGRE TB 880 Table 2-3)
rho_filler = 6.0                       -- Fillers  (CIGRE TB 880 Table 12-11)
rho_core_binder = 6.0                  -- Binding tape (CIGRE TB 880 Note 2 - Table 12-11)
rho_armour = 1.0/45                    -- Armour (assumed)
rho_outer_covering = 6.0               -- Outer covering (CIGRE TB 880 Table 2-3)
rho_seabed = seabed_thermal_resistivity

-- Calculated radius and thicknesses
r_conductor = conductor_diameter / 2
r_conductor_tape = conductor_tape_diameter / 2
r_conductor_screen = conductor_screen_outer_diameter / 2
r_insulation = insulation_outer_diameter / 2
r_insulation_screen = insulation_screen_outer_diameter / 2
r_water_blocking = water_blocking_outer_diameter / 2
r_lead_sheath = lead_sheath_outer_diameter / 2
r_pe_sheath = pe_sheath_outer_diameter / 2
r_filler = filler_diameter / 2
r_cable_core_binder = cable_core_binder_diameter / 2
r_armour = armour_outer_diameter / 2
r_cable_outer = cable_outer_diameter / 2

-- New Heat Flow document, units in millimeters, planar, default precision
newdocument(2)  -- 2 = Heat Flow problem in FEMM
depth_mm = 1 -- depth = 1 mm for per-unit-depth results
min_angle = 30 -- Related to meshing
hi_probdef("millimeters", "planar", 1e-8, depth_mm, min_angle)  -- depth = 1 mm for per-unit-depth results
path = "C:\\Desktop\\Thermal Simulations\\" --Change the path here to store the file
showconsole()

-- -------------------------
-- Loss calculations
-- -------------------------
-- Calculate layup factor for three-core stranded arrangement
-- Rationale: Accounts for increased path length due to stranding (IEC 60287 assumes straight, TB 880 GP23 corrects).
function calculate_core_layup_factor()
  local factor_Cfl = 1.29 -- Based on IEC 60287 Table A.1 for 3 number of cores
  local layup_term = pi * (factor_Cfl * pe_sheath_outer_diameter) / core_lay_length
  local f_layup = sqrt(1 + layup_term^2)
    
  print("\n=== CORE LAYUP EFFECT CALCULATION ===")
  print("Layup factor: " .. f_layup)
  return f_layup
end

-- Conductor AC resistance calculation 
function calculate_conductor_losses(conductor_temp)
  local factor_due_to_armour = "None"
  -- Factor related to non-magnetic and magnetic armour wires - CIGRE TB 880 Guidance Point 38 for magnetic armour wires and IEC 60287-1-1 Section 2.1.1
  if material_armour == "Galvanized Stainless Steel" then
    factor_due_to_armour = 1
  elseif material_armour == "Galvanized Steel" then
    factor_due_to_armour = 1.5
  else
    print("Armour material is wrongly set")
  end

  -- DC resistance at operating temperature - IEC 60287 1-1 section 2.1.1
  local R_dc = R_dc_20 * (1 + alpha * (conductor_temp - 20))
  
  -- Skin effect calculation - IEC 60287-1-1 Section 2.1.2
  local x_s_squared = (8 * pi * system_frequency * 1e-7 / R_dc) * skin_effect_factor
  local x_s = sqrt(x_s_squared)
  
  local y_s = "None"
  if x_s <= 2.8 then
    y_s = x_s^4 / (192 + 0.8 * x_s^4)
  elseif x_s <= 3.8 then
    y_s = -0.136 - 0.0177 * x_s + 0.0563 * x_s * x_s
  else
    y_s = 0.354 * x_s - 0.733 
  end
  
  -- Proximity effect calculation - IEC 60287-1-1 Section 2.1.4
  local x_p_squared = (8 * pi * system_frequency * 1e-7 / R_dc) * proximity_effect_factor
  local x_p = sqrt(x_p_squared)
  local dc_over_s = conductor_diameter / phase_spacing
  local y_p = (x_p^4 / (192 + 0.8 * x_p^4)) * (dc_over_s^2) * (0.312 * (dc_over_s^2) + (1.18 / ((x_p^4 / (192 + 0.8 * x_p^4)) + 0.27))) 

  -- AC resistance per meter of core - IEC 60287-1-1 Section 2.1
  local R_ac = R_dc * (1 + factor_due_to_armour *( y_s + y_p))
  
  -- Adjust for layup factor in case in need to be added later
  --R_ac = layup_factor * R_ac -- Attention here 
      
  return R_ac
end

-- Lead sheath losses - IEC 60287 Section 5.3
function calculate_sheath_losses(sheath_temp , R_con_ac)
  -- Lead sheath thickness calculation
  local lead_sheath_thickness = (lead_sheath_outer_diameter - water_blocking_outer_diameter) / 2

  local factor_due_to_armour = "None"
  -- Factor related to non-magnetic and magnetic armour wires - IEC 60287-1-1 Section 5.3.2 for non-magnetic armour and CIGRE TB 880 Guidance Point 38 for magnetic armour wires
  if material_armour == "Galvanized Stainless Steel" then
    factor_due_to_armour = 1
  elseif material_armour == "Galvanized Steel" then
    factor_due_to_armour = 1.5
  else
    print("Armour material is wrongly set")
  end

  -- Sheath cross-sectional area
  local A_s = pi * lead_sheath_thickness * (lead_sheath_outer_diameter - lead_sheath_thickness)  -- mm²

  -- Sheath resistance at operating temperature - IEC 60287-1-1 Section 5.3.1
  local resistivity_sheath = resistivity_sheath_20 * (1 + alpha_sheath * (sheath_temp - 20))
  local R_s = (resistivity_sheath * layup_factor) / (A_s * 1e-6)  -- Ω/m 

  -- Calculation of the reactance X per unit length of the sheath - IEC 60287 Section 5.3.2
  local mean_diam_sheath = lead_sheath_outer_diameter - lead_sheath_thickness  
  local omega = 2 * pi * system_frequency
  local X_core =  2 * omega * 1e-7 * log(2 * phase_spacing / mean_diam_sheath ) * layup_factor

  -- Calculation of lambda 1' for circulating losses of the sheath - IEC 60287-1-1 Section 5.3.2 for non-magnetic armour and CIGRE TB 880 Guidance Point 38 for magnetic armour wires
  local RS_X = R_s / X_core 
  local lambda_1_circulating = factor_due_to_armour * R_s / ( R_con_ac * (1+ RS_X * RS_X )) 

  -- Calculation of lambda 1" for eddy currents of the sheath - IEC 60287-1-1 Section 5.3.2  - This approach differs from IEC 60287-1-1 section 5.3.2. Refer to CIGRE TB 880 Guidance Point 6
  local m = omega * 1e-7 / R_s
  local lambda_0 = 3.0 * ((m*m) / (1.0 + m*m)) * ((mean_diam_sheath / (2.0*phase_spacing))^2)
  local delta_1 = (1.14 * (m^2.45) + 0.33) * ((mean_diam_sheath / (2.0 * phase_spacing)) ^ (0.92 * m + 1.66))
  local delta_2 = 0
  local b1 = sqrt( 4 * pi * omega * 1e-7 / resistivity_sheath)
  local gs = 1 + (lead_sheath_thickness/lead_sheath_outer_diameter)^1.74 * (b1 * lead_sheath_outer_diameter * 1e-3 - 1.6)
  local lambda_1_eddy = (R_s / R_con_ac) * ( gs * lambda_0 * (1.0 + delta_1 + delta_2) + ((b1 * lead_sheath_thickness)^4) / (12.0 * 1e12) )

  -- Effect of eddy currents Cf-factor - CIGRE TB 880 Guidance Point 31 - IEC 60287-1-1 Section 5.3.6
  local M = RS_X
  local N = RS_X
  local F = (4.0 * M * M * N * N + (M + N)^2) / (4.0 * (M * M + 1.0) * (N * N + 1.0))

  -- Calculation of total loss factor - IEC 60287-1-1 Section 5.3.1
  local lambda_1 = lambda_1_circulating + F * lambda_1_eddy
  
  return lambda_1, R_s, lambda_1_circulating
end

-- Calculation of Armour losses  - CIGRE TB 880 Section 
function calculate_armour_losses(armour_temp, R_con_ac, R_sheath, lambda_1_circ)
  local lambda_2 = "None"
  local mean_diam_armour = armour_outer_diameter - diameter_armour_wires
  -- Angular frequency  
  local omega = 2 * pi * system_frequency

  if material_armour == "Galvanized Stainless Steel" then
    lambda_2 = 0 -- CIGRE TB 880 Guidance Point 36 
  elseif material_armour == "Galvanized Steel" then
    -- Calculation of cross-sectional area of the armour
    local A_armour = number_armour_wires * pi * (diameter_armour_wires/2)^2 --mm²
    -- Ratio of the lenght of the wires to the length of the cable
    local layup_armour = sqrt(1+(pi * mean_diam_armour / armour_lay_length)^2)
    -- Resistance of Armour at 20oC 
    local R_AO_temp = resistivity_armour_20 / (A_armour * 1e-6)  -- Convert A_armour mm² to m²
    local k = (1.4-1.2)/(5-2) * (diameter_armour_wires - 2) + 1.2-- Guidance Point 34 and IEC 60287-1-1 Section 5.4.3.1
    local R_AO = k * R_AO_temp * layup_armour
    local R_A = R_AO * (1 + alpha_armour * (armour_temp - 20))
    -- Calculation of the value c - CIGRE TB 880 Guidance Point 33
    local C = filler_diameter / 2 - pe_sheath_outer_diameter / 2 --mm²
    -- Calculation of factor lamba_2 of the armour - CIGRE TB 880 Guidance Point 37 - IEC 60287-1-1 Section 5.4.3.5
    local factor_lamba2 = (1 - R_con_ac * lambda_1_circ / R_sheath)
    if factor_lamba2 < 0 then 
        factor_lamba2 = 0
    end
    -- Calculation of lamba_2
    lambda_2 = 1.23 * (R_A / R_con_ac) * ((2 * C * 1e-3) / (mean_diam_armour * 1e-3))^2 * (1 / (((2.77 * R_A * 1e6) / omega)^2 + 1)) * factor_lamba2 -- IEC 60287-1-1 Section 5.4.3.3.1
  else
    print("Armour material is wrongly set")
  end

  return lambda_2
end

-- -------------------------
-- Helpers
-- -------------------------
function con_from_rho(rho_Km_per_W)
  -- Convert rho [K·m/W] to k [W/(mm·K)] with k = (1/rho) / 1000
  local k_W_per_mK = 1.0 / rho_Km_per_W
  return k_W_per_mK 
end


function add_full_circle(xc, yc, r, maxsegdeg)
  -- Build a full circle from two 180° arcs between (xc±r,yc)
  hi_addnode(xc + r, yc)
  hi_addnode(xc - r, yc)
  -- Top semicircle: right -> left, 180 degrees
  hi_addarc(xc + r, yc, xc - r, yc, 180.0, maxsegdeg)
  -- Bottom semicircle: left -> right, 180 degrees
  hi_addarc(xc - r, yc, xc + r, yc, 180.0, maxsegdeg)
end

function add_annulus_label(xc, yc, r_in, r_out)
  -- Place a label radially at the midpoint, on +x axis
  local r_mid = 0.5 * (r_in + r_out)
  hi_addblocklabel(xc + r_mid, yc)
end

function set_conductor_on_circle(xc, yc, r, conductor_name, maxsegdeg)
  hi_selectarcsegment(xc, yc + r)
  hi_setarcsegmentprop(maxsegdeg, "<None>", 0, 0, conductor_name)
  hi_clearselected()
  hi_selectarcsegment(xc, yc - r)
  hi_setarcsegmentprop(maxsegdeg, "<None>", 0, 0, conductor_name)
  hi_clearselected()
end
-- -------------------------
-- Define materials once (isotropic kx=ky, zero qv, zero volumetric heat capacity)
-- -------------------------
function define_cable_materials()
  -- Compute k in W/(mm·K)
  local k_conductor          = con_from_rho(rho_conductor)
  local k_conductor_tape     = con_from_rho(rho_conductor_tape)
  local k_conductor_screen   = con_from_rho(rho_conductor_screen)
  local k_insulation         = con_from_rho(rho_xlpe_insulation)
  local k_insulation_screen  = con_from_rho(rho_insulation_screen)
  local k_water_blocking     = con_from_rho(rho_water_blocking)
  local k_lead_sheath        = con_from_rho(rho_lead_sheath)
  local k_pe                 = con_from_rho(rho_pe_sheath)
  local k_filler             = con_from_rho(rho_filler)
  local k_core_binder        = con_from_rho(rho_core_binder)
  local k_armour             = con_from_rho(rho_armour)
  local k_cable_outer        = con_from_rho(rho_outer_covering)
  local k_seabed             = con_from_rho(rho_seabed)  

  -- Add FEMM heat-flow materials: hi_addmaterial(name, kx, ky, qv, kt)
  -- qv = 0 (no volumetric generation by default), kt = 0 (no capacity in steady state)
  hi_addmaterial("Conductor",            k_conductor,        k_conductor,        0, 0)
  hi_addmaterial("Conductor Tape",       k_conductor_tape,   k_conductor_tape,   0, 0)
  hi_addmaterial("Conductor Screen",     k_conductor_screen, k_conductor_screen, 0, 0)
  hi_addmaterial("XLPE Insulation",      k_insulation,       k_insulation,       0, 0)
  hi_addmaterial("Insulation Screen",    k_insulation_screen,k_insulation_screen,0, 0)
  hi_addmaterial("Water Blocking Tape",  k_water_blocking,   k_water_blocking,   0, 0)
  hi_addmaterial("Lead Sheath",          k_lead_sheath,      k_lead_sheath,      0, 0)
  hi_addmaterial("PE Sheath",            k_pe,               k_pe,               0, 0)
  hi_addmaterial("Fillers",              k_filler,           k_filler,           0, 0)
  hi_addmaterial("Cable core binder",    k_core_binder,      k_core_binder,      0, 0)
  hi_addmaterial("Armour",               k_armour,           k_armour,           0, 0)
  hi_addmaterial("PP Outer Layer",       k_cable_outer,      k_cable_outer,      0, 0)
  hi_addmaterial("Seabed Soil",          k_seabed,           k_seabed,           0, 0)
end

----------------------------
-- Main function requested: draw and assign properties around (xc, yc)
function draw_power_core_cable(xc, yc, maxsegdeg, cable_id, core_id)

  -- Interfaces outwards
  add_full_circle(xc, yc, r_conductor, maxsegdeg)
  add_full_circle(xc, yc, r_conductor_tape, maxsegdeg)
  add_full_circle(xc, yc, r_conductor_screen, maxsegdeg)
  add_full_circle(xc, yc, r_insulation, maxsegdeg)
  add_full_circle(xc, yc, r_insulation_screen, maxsegdeg)
  add_full_circle(xc, yc, r_water_blocking, maxsegdeg)
  add_full_circle(xc, yc, r_lead_sheath, maxsegdeg)
  add_full_circle(xc, yc, r_pe_sheath, maxsegdeg)

  -- Place labels for each region (centered radially within each annulus)
  -- Conductor solid
  local cond_name = "Conductor" .. cable_id .."Core" ..core_id
  hi_addconductorprop(cond_name, 0, conductor_qc_per_core, 0)  -- Type 0: prescribed heat flow
  set_conductor_on_circle(xc, yc, r_conductor, cond_name, maxsegdeg)
  hi_addblocklabel(xc, yc)
  hi_selectlabel(xc, yc)
  hi_setblockprop("Conductor", 1, 0, 0)
  hi_clearselected()

  -- Conductor Tape: 
  add_annulus_label(xc, yc, r_conductor, r_conductor_tape)
  hi_selectlabel(xc + 0.5*(r_conductor + r_conductor_tape), yc)
  --local meshsize_con_tap = (r_conductor_tape - r_conductor) / 3 -- To solve Triangle issues by setting meshsize in thin layers
  hi_setblockprop("Conductor Tape", 1, 0, 0)
  hi_clearselected()

  -- Conductor Screen:
  add_annulus_label(xc, yc, r_conductor_tape, r_conductor_screen)
  hi_selectlabel(xc + 0.5*(r_conductor_tape + r_conductor_screen), yc)
  --local meshsize_con_sc = (r_conductor_screen - r_conductor_tape) / 3 -- To solve Triangle issues by setting meshsize in thin layers
  hi_setblockprop("Conductor Screen", 1, 0, 0)
  hi_clearselected()

  -- XLPE Insulation: 
  add_annulus_label(xc, yc, r_conductor_screen, r_insulation)
  hi_selectlabel(xc + 0.5*(r_conductor_screen + r_insulation), yc)
  hi_setblockprop("XLPE Insulation", 1, 0, 0)
  hi_clearselected()

  -- Insulation Screen:
  add_annulus_label(xc, yc, r_insulation, r_insulation_screen)
  hi_selectlabel(xc + 0.5*(r_insulation + r_insulation_screen), yc)
  --local meshsize_ins_sc = (r_insulation_screen - r_insulation) / 3 -- To solve Triangle issues by setting meshsize in thin layers
  hi_setblockprop("Insulation Screen", 1, 0, 0)
  hi_clearselected()

  -- Water Blocking Tape: 
  add_annulus_label(xc, yc, r_insulation_screen, r_water_blocking)
  hi_selectlabel(xc + 0.5*(r_insulation_screen + r_water_blocking), yc)
  --local meshsize_water_tape = (r_water_blocking - r_insulation_screen) / 3 -- To solve Triangle issues by setting meshsize in thin layers
  hi_setblockprop("Water Blocking Tape", 1, 0, 0)
  hi_clearselected()

  -- Lead Sheath: 
  local sheath_name = "Sheath" .. cable_id.. "Core" ..core_id
  hi_addconductorprop(sheath_name, 0, sheath_qc_per_core, 0)
  set_conductor_on_circle(xc, yc, r_water_blocking, sheath_name, maxsegdeg)  -- Inner
  set_conductor_on_circle(xc, yc, r_lead_sheath, sheath_name, maxsegdeg)  -- Outer
  add_annulus_label(xc, yc, r_water_blocking, r_lead_sheath)
  hi_selectlabel(xc + 0.5*(r_water_blocking + r_lead_sheath), yc)
  hi_setblockprop("Lead Sheath", 1, 0, 0)
  hi_clearselected()

  -- PE Sheath: 
  add_annulus_label(xc, yc, r_lead_sheath, r_pe_sheath)
  hi_selectlabel(xc + 0.5*(r_lead_sheath + r_pe_sheath), yc)
  hi_setblockprop("PE Sheath", 1, 0, 0)
  hi_clearselected()

end

function draw_submarine_cable (x_center, y_center,cable_id )

  -- Draw all circular interfaces; fine arc segmentation for good annulus meshes
  local maxsegdeg = 1  -- degrees per element on arcs

  --Cable Cores
  -- Angles for trefoil: 0° (top), 120° (bottom-left), 240° (bottom-right)
  local angles = {0, 120, 240}
  -- For an equilateral triangle of side = outer_diam, circumradius = side / sqrt(3)
  local R = pe_sheath_outer_diameter / sqrt(3)

  -- To solve the meshing issue
  local vertical_shift = 0.01 --in mm
  local horizontal_shift = 0.01 --in mm

  -- Place three centers and draw the three cables
  -- Mapping 0° -> top, 120° -> bottom-left, 240° -> bottom-right:
  -- use dx = -R * sin(a), dy = R * cos(a) to align 0° with "up"
  for i = 1, 3 do
    local a = angles[i]
    local a_rad = a * (pi / 180) 
    local dx = -R * sin(a_rad)
    local dy =  R * cos(a_rad)
    local xc = x_center + dx
    local yc = 0
    if i == 1 then
      yc = y_center + dy + vertical_shift
    else
      yc = y_center + dy - vertical_shift
    end
    if i == 3 then  -- Apply to one core, e.g., bottom-left
      xc = xc + horizontal_shift
    elseif i == 2 then  -- Apply opposite to bottom-right for balance
      xc = xc - horizontal_shift
    end
    draw_power_core_cable(xc, yc, maxsegdeg,cable_id, i)
  end

  -- After Cable cores
  add_full_circle(x_center, y_center, r_filler, maxsegdeg) 
  add_full_circle(x_center, y_center, r_cable_core_binder, maxsegdeg)
  add_full_circle(x_center, y_center, r_armour, maxsegdeg)
  add_full_circle(x_center, y_center, r_cable_outer, maxsegdeg)

  -- Fillers: 
  add_annulus_label(x_center, y_center, 0, 0)
  hi_selectlabel(x_center, y_center)
  hi_setblockprop("Fillers", 1, 0, 0)
  hi_clearselected()

  -- Cable Binder Core:  To solve the problem with meshing
  add_annulus_label(x_center, y_center, r_filler, r_cable_core_binder)
  hi_selectlabel(x_center + 0.5*(r_filler + r_cable_core_binder), y_center)
  --local meshsize_core_bin = (r_cable_core_binder - r_filler) / 3 -- To solve Triangle issues by setting meshsize in thin layers
  hi_setblockprop("Cable core binder", 1, 0, 0)
  hi_clearselected()

  -- Armour: 
  local armour_name = "Armour"..cable_id
  hi_addconductorprop(armour_name, 0, armour_qc, 0)
  set_conductor_on_circle(x_center, y_center, r_cable_core_binder, armour_name, maxsegdeg)  -- Inner
  set_conductor_on_circle(x_center, y_center, r_armour, armour_name, maxsegdeg)  -- Outer
  add_annulus_label(x_center, y_center, r_cable_core_binder, r_armour)
  hi_selectlabel(x_center + 0.5*(r_cable_core_binder + r_armour), y_center)
  hi_setblockprop("Armour", 1, 0, 0)
  hi_clearselected()

  -- Cable Outer Layer: 
  add_annulus_label(x_center, y_center, r_armour, r_cable_outer)
  hi_selectlabel(x_center + 0.5*(r_armour + r_cable_outer), y_center)
  hi_setblockprop("PP Outer Layer", 1, 0, 0)
  hi_clearselected()
end

function domain(x_left, x_right, y_bottom, y_top )
  -- Draw rectangular soil domain
  hi_addnode(x_left, y_top)
  hi_addnode(x_right, y_top)
  hi_addnode(x_right, y_bottom)
  hi_addnode(x_left, y_bottom)
  hi_addsegment(x_left, y_top, x_right, y_top)
  hi_addsegment(x_right, y_top, x_right, y_bottom)
  hi_addsegment(x_right, y_bottom, x_left, y_bottom)
  hi_addsegment(x_left, y_bottom, x_left, y_top)

  -- Define boundary property for seabed surface (fixed temperature)
  hi_addboundprop("Seabed Surface", 0, seabed_temperature+273, 0, 0, 0)  -- Type 0: Fixed T. FEMM use Kelvin so 273 needs to be added

  -- Apply BC to top segment
  local mid_x = (x_left + x_right) / 2
  hi_selectsegment(mid_x, y_top)
  hi_setsegmentprop("Seabed Surface", 0, 0, 0, 0, 0)
  hi_clearselected()

  -- Place soil block label 
  local label_x = 0
  local mean_y = 0.5 * (burial_depth - r_cable_outer)  
  local label_y = -burial_depth + mean_y
  hi_addblocklabel(label_x, label_y)
  hi_selectlabel(label_x, label_y)
  hi_setblockprop("Seabed Soil", 1, 0, 0)  -- Automesh enabled
  hi_clearselected()
end

-- -------------------------
-- Problem initialization (run once per model)
-- -------------------------
layup_factor = calculate_core_layup_factor()

-- Accounting for lay lenght factor
rho_conductor_tape = rho_conductor_tape / layup_factor
rho_conductor_screen = rho_conductor_screen / layup_factor
rho_xlpe_insulation = rho_xlpe_insulation / layup_factor
rho_insulation_screen = rho_insulation_screen / layup_factor
rho_water_blocking = rho_water_blocking / layup_factor
rho_pe_sheath = rho_pe_sheath / layup_factor

-- Set of temperatures as per ambient temperature
Conductor_temp = seabed_temperature
Sheath_temp = seabed_temperature
Armour_temp = seabed_temperature

-- Initial Losses 
qc_scale = depth_mm / 1000  -- W/m to W per mm depth (depth of problem and convertion from m to mm)
qc_sign = 1.0 
R_con_ac = calculate_conductor_losses(Conductor_temp)
lambda_1, R_sheath, lambda_1_circulating = calculate_sheath_losses(Sheath_temp , R_con_ac)
lambda_2 = calculate_armour_losses(Armour_temp, R_con_ac, R_sheath, lambda_1_circulating)
Conductor_losses = Current^2 * R_con_ac
Sheath_losses = Current^2 * R_con_ac * lambda_1
Armour_losses = Current^2 * R_con_ac * lambda_2 * 3 -- Armour losses should account all cables
conductor_qc_per_core = qc_sign * Conductor_losses * qc_scale
sheath_qc_per_core = qc_sign * Sheath_losses * qc_scale
armour_qc = qc_sign * Armour_losses * qc_scale

-- Burial positioning
y_top = 0           -- Seabed surface (mm)
y_center = y_top-burial_depth  -- Shift cables down to burial depth (mm)

-- Define soil domain dimensions (adjust as needed for convergence)
domain_width_margin = 10000  -- margin on each side (mm)
domain_depth_below = 20000   -- below lowest cable point (mm)
y_bottom = y_center - r_cable_outer - domain_depth_below -- Seabed surface (mm)

define_cable_materials()

if number_cables == 1 then
  draw_submarine_cable (0,y_center,1)
  x_left = -r_cable_outer - domain_width_margin 
  x_right = r_cable_outer + domain_width_margin
elseif number_cables == 2 then
  draw_submarine_cable (0,y_center,1)
  draw_submarine_cable (axial_separation,y_center,2)
  x_left = -r_cable_outer - domain_width_margin 
  x_right = r_cable_outer + domain_width_margin + axial_separation 
elseif number_cables == 3 then
  draw_submarine_cable (0,y_center,1)
  draw_submarine_cable (axial_separation,y_center,2)
  draw_submarine_cable (-axial_separation,y_center,3)
  x_left = -r_cable_outer - domain_width_margin - axial_separation 
  x_right = r_cable_outer + domain_width_margin + axial_separation 
else
  print("Number of cables error")
end

domain(x_left, x_right, y_bottom, y_top )

-- Save the initial problem
hi_saveas(path.."thermal_calc.feh") -- Change the path here

-- Iteration loop for convergence
tolerance = 0.001
max_iter = 50
converged = false
iter = 0
num_cores = number_cables * 3 -- 3 core cables 
num_armours = number_cables

while not converged and iter < max_iter do
  iter = iter + 1
  print("Iteration " .. iter .. ": Conductor Temp = " .. Conductor_temp .. ", Sheath Temp = " .. Sheath_temp .. ", Armour Temp = " .. Armour_temp)

  -- Update conductor properties with current qc values
  for cable_id = 1, number_cables do
  for core_id = 1, 3 do
    local cond_name = "Conductor" .. cable_id .. "Core" .. core_id
    hi_modifyconductorprop(cond_name, 2, conductor_qc_per_core)
    local sheath_name = "Sheath" .. cable_id .. "Core" .. core_id
    hi_modifyconductorprop(sheath_name, 2, sheath_qc_per_core)
  end
  local armour_name = "Armour" .. cable_id
  hi_modifyconductorprop(armour_name, 2, armour_qc)
  end

  -- Analyze and load solution
  hi_analyze(1)
  hi_loadsolution()

  -- Extract average temperatures from postprocessor using conductor properties
  local sum_cond_t = 0
  for cable_id = 1, number_cables do
    for core_id = 1, 3 do
      local cond_name = "Conductor" .. cable_id .. "Core" .. core_id
      local temp, flux = ho_getconductorproperties(cond_name)
      sum_cond_t = sum_cond_t + ( temp - 273 ) -- convection from Kelvin to Celsius
    end
  end
  local new_cond_t = sum_cond_t / num_cores

  local sum_sheath_t = 0
  for cable_id = 1, number_cables do
    for core_id = 1, 3 do
      local sheath_name = "Sheath" .. cable_id .. "Core" .. core_id
      local temp, flux = ho_getconductorproperties(sheath_name)
      sum_sheath_t = sum_sheath_t + ( temp - 273 ) -- convection from Kelvin to Celsius
    end
  end
  local new_sheath_t = sum_sheath_t / num_cores

  local sum_arm_t = 0
  for cable_id = 1, number_cables do
    local armour_name = "Armour" .. cable_id
    local temp, flux = ho_getconductorproperties(armour_name)
    sum_arm_t = sum_arm_t + ( temp - 273 ) -- convection from Kelvin to Celsius
  end
  local new_arm_t = sum_arm_t / num_armours

  -- Check convergence
  converged = abs(new_cond_t - Conductor_temp) < tolerance and abs(new_sheath_t - Sheath_temp) < tolerance and abs(new_arm_t - Armour_temp) < tolerance

  -- Update temperatures
  Conductor_temp = new_cond_t
  Sheath_temp = new_sheath_t
  Armour_temp = new_arm_t

  -- Recalculate losses and qc based on new temperatures
  R_con_ac = calculate_conductor_losses(Conductor_temp)
  lambda_1, R_sheath, lambda_1_circulating = calculate_sheath_losses(Sheath_temp , R_con_ac)
  lambda_2 = calculate_armour_losses(Armour_temp, R_con_ac, R_sheath, lambda_1_circulating)
  Conductor_losses = Current^2 * R_con_ac
  Sheath_losses = Current^2 * R_con_ac * lambda_1
  Armour_losses = Current^2 * R_con_ac * lambda_2 * 3 -- Armour losses should account all cables
  conductor_qc_per_core = qc_sign * Conductor_losses * qc_scale
  sheath_qc_per_core = qc_sign * Sheath_losses * qc_scale
  armour_qc = qc_sign * Armour_losses * qc_scale
end

if converged then
  print("Converged after " .. iter .. " iterations.")
else
  print("Did not converge within " .. max_iter .. " iterations.")
end

print("Final Conductor Temperature: " .. Conductor_temp .. " °C")
print("Final Conductor Losses: " .. Conductor_losses .. " W/m per core")
print("Final Sheath Losses: " .. Sheath_losses .. " W/m per core")
print("Final Armour Losses: " .. Armour_losses .. " W/m per cable")


