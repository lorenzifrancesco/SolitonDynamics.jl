function sim_info(sim)
  @unpack_Sim N, L, manual = sim
  print("\nSimulation in $(length(N))D, with a number of steps $(N)\n")
end

"""
Integration check
"""
function check_module()
  print("Module is calling correclty")
end