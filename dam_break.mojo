
from utils.static_tuple import StaticTuple
from utils.vector import InlinedFixedVector
from python import Python

# Meshing and Geometry
alias lenght: Float64 = 1
alias height: Float64 = 1
alias mesh_size: Int = 64

# Solution proceedure
var cfl: Float64 = 0.5
var time_step: Int = 0
let time_end: Float64 = 0
var time_n: Float64 = 0

# Field Data
var volume_fraction: InlinedFixedVector[mesh_size*mesh_size, Float64] = 0
var pressure: InlinedFixedVector[mesh_size*mesh_size, Float64] = 0
var velocity: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0

var gradient_volume_fraction: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0
var gradient_pressure: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, Float64]] = 0
var gradient_velocity: InlinedFixedVector[mesh_size*mesh_size, StaticTuple[2, StaticTuple[2, Float64]]] = 0

fn make_mesh():
  return

fn compute_stability() -> Float64:

  return 0

fn compute_gradients():
  return 

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
  make_mesh()
  print("solving...")
  while(not solve_time_step()): # so it does not inf loop - need stabiltiy
    print("time step ",  time_n)  

  volume_fraction[0] = 5
  
  print(volume_fraction[0])
