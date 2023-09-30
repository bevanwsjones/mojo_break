
from utils.static_tuple import StaticTuple
from utils.vector import InlinedFixedVector
import math
from python import Python



"""
  Grid indexing:

  m - number of cells in a cardinal sirection (m)
  c_m - number of cells (m*m)
  f_m - maximum number of faces (m + 1)*(m + 1)

                            f = (m*(m + 1) + m + i)
                        |---------|
                        |         |
              f = i - 1 |    i    | f = i 
                        |         |
                        |---------|
                            f = (m*(m + 1) + i)

  first horizontal face -> f_h0 = m*(m + 1)

boundary face key
0 - not a boundary face
1 - left boundary
2 - right boundary
4 - bottom boundary
5 - top boundary

"""


# Meshing and Geometry
alias lenght: Float64 = 1
alias height: Float64 = 1
alias mesh_size: Int = 64

let face_h0 : Int = mesh_size*(mesh_size + 1)
let face_m : Int = (mesh_size + 1)*(mesh_size + 1)

# Solution proceedure
var cfl: Float64 = 0.5
var time_step: Int = 0
let time_end: Float64 = 0
var time_n: Float64 = 0

# Geometric Data
let delta_x : Float64 = lenght/Float64(mesh_size)
let delta_y : Float64 = height/Float64(mesh_size)
let face_area_x : Float64 = height/Float64(mesh_size)
let face_area_y : Float64 = lenght/Float64(mesh_size)
let volume = 1.0/(delta_x*delta_y)
var boundary_face_type : InlinedFixedVector[mesh_size*mesh_size, Int8] = 0


# material properties
let fluid_density : StaticTuple[2, Float64] = VariadicList[Float64](1000.0, 1.0)
let fluid_viscosity : StaticTuple[2, Float64] = VariadicList[Float64](1.0, 1.0e-3)

# Field Data
let N : Int = 1

var volume_fraction: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0
var pressure: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0
var velocity: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, StaticTuple[2, Float64]]] = 0

var gradient_volume_fraction: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0
var gradient_pressure: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0
var gradient_velocity: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, StaticTuple[2, Float64]]] = 0

# ----------------------------------------------------------------------------------------------------------------------
# Geometic Functions
# ----------------------------------------------------------------------------------------------------------------------

fn cell_coordinate(i_x : Int, i_y : Int) -> StaticTuple[2, Float64]:
  """
  
  """

  return VariadicList(delta_x*(0.5 + Float64(i_x)), delta_y*(0.5 + Float64(i_y)))

fn cell_0(i_face : Int) -> Int:
  return (i_face - 1) if i_face < face_h0 else (mesh_size*(mesh_size + 1) + i_face) 

fn cell_1(i_face : Int) -> Int:
  return (i_face - 1) if i_face < face_h0 else (mesh_size*(mesh_size + 1) + mesh_size + i_face) 

fn is_cell_0_first(i_face : Int) -> Bool:
  return i_face < face_h0 and (i_face%mesh_size != 0)

fn is_cell_0_last(i_face : Int) -> Bool:
  return i_face > (face_h0 + mesh_size)

fn face_area_normal(i_face : Int) -> StaticTuple[2, Float64]:
  return VariadicList[Float64](face_area_y, 0.0) if i_face < face_h0 else VariadicList[Float64](0.0, face_area_x) 

fn is_bounday_face(i_face : Int) -> Int8:
  if i_face % (mesh_size + 1) == 0: return 1
  elif i_face % mesh_size == 0: return 2
  elif face_h0 <= i_face < (face_h0 + mesh_size): return 3 
  elif face_m - mesh_size < i_face: return 4
  else: return 0

# ----------------------------------------------------------------------------------------------------------------------
# Support Functions and discretisation
# ----------------------------------------------------------------------------------------------------------------------

fn unit_lerp(borrowed x_i : Float64, borrowed y_0 : Float64, borrowed y_1 : Float64) -> Float64:
  return y_0*x_i + (1.0 - x_i)*y_1

fn density(borrowed volume_fraction : Float64) -> Float64:
  let clamped_volume_faction : Float64 = math.clamp(volume_fraction, 0.0, 1.0)
  return unit_lerp(clamped_volume_faction, fluid_density[0], fluid_density[1])

fn viscosity(borrowed volume_fraction : Float64) -> Float64:
  let clamped_volume_faction : Float64 = math.clamp(volume_fraction, 0.0, 1.0)
  return unit_lerp(clamped_volume_faction, fluid_density[0], fluid_density[1])

fn arithmetic_mean(value_0 : Float64, value_1 : Float64) -> Float64:
  return 0.5*(value_0 + value_1)

# ----------------------------------------------------------------------------------------------------------------------
# Stability
# ----------------------------------------------------------------------------------------------------------------------

fn compute_stability() -> Float64:

  return 0


# ----------------------------------------------------------------------------------------------------------------------
# Stability
# ----------------------------------------------------------------------------------------------------------------------

fn compute_gradients():
  
  for i_face in range(face_m):
    let i_cell_0 = cell_0(i_face)
    let i_cell_1 = cell_1(i_face)

    if is_cell_0_first(i_face):
      gradient_volume_fraction[i_cell_0][0] *= 0.0
      gradient_volume_fraction[i_cell_0][1] *= 0.0

    let area_normal = face_area_normal(i_face)
    if boundary_face_type[i_face]:
      let i_cell = i_cell_0 if i_face % 2 != 0 else i_cell_1
      let coef_sign = -1.0 if i_face % 2 != 0 else 1.0
      gradient_volume_fraction[i_cell][0] += coef_sign*volume_fraction[N][i_cell]*area_normal[0]
      gradient_volume_fraction[i_cell][1] += coef_sign*volume_fraction[N][i_cell]*area_normal[1]
    else:
      let face_value = arithmetic_mean(volume_fraction[N][i_cell_0], volume_fraction[N][i_cell_1])
      gradient_volume_fraction[i_cell_0][0] -= face_value*area_normal[0]
      gradient_volume_fraction[i_cell_0][1] -= face_value*area_normal[1]
      gradient_volume_fraction[i_cell_1][0] += face_value*area_normal[0]
      gradient_volume_fraction[i_cell_1][1] += face_value*area_normal[1]

    if is_cell_0_last(i_face):
      gradient_volume_fraction[i_cell_0][0] /= volume
      gradient_volume_fraction[i_cell_0][1] /= volume

fn compute_gradients_pressure():
  
  return 

# ----------------------------------------------------------------------------------------------------------------------
# Equation Solver Function
# ----------------------------------------------------------------------------------------------------------------------

fn step_vof(delta_time: Float64):
  return

fn step_0(delta_time: Float64):
  return

fn step_1(delta_time: Float64):
  
  return

fn step_2(delta_time: Float64):
  return

fn solve_time_step() -> Bool:

  var delta_time = compute_stability()
  time_n += delta_time

  compute_gradients()
  
  step_vof(time_step)

  step_0(time_step)

  step_1(time_step)

  step_2(time_step)

  time_step += 1 
  return time_n >= time_end

fn main():

  for i_face in range(face_m):
    boundary_face_type[i_face] = is_bounday_face(i_face)

  print("solving...")
  while(not solve_time_step()): # so it does not inf loop - need stabiltiy
    print("time step ",  time_n)  

