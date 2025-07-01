rotateMapCorrectly = function(map, optimisations_to_rotate = 1){
  
  getMapAngle = function(map, ag1, ag2, op){
    coords = Racmacs::agCoords(map, optimization_number = op)[c(ag1, ag2),]
    x = atan(diff(coords[,2]) / diff(coords[,1])) * 180 / pi
    if (coords[2,1]> coords[1,1]) x = 180 + x
    x
  }
  
  for (op in optimisations_to_rotate){
    
    map_rot = Racmacs::rotateMap(
      map,
      -(180-getMapAngle(map, "Wu1", "XBB.1.5", op = op)),
      optimization_number = op
    )
    
    if (agCoords(map_rot)["Wu1",2] < agCoords(map_rot)["JN.1",2]){
      map_ref = Racmacs::reflectMap(map, optimization_number = op)
      map_rot = Racmacs::rotateMap(
        map_ref,
        -(180-getMapAngle(map_ref, "Wu1", "XBB.1.5", op = op)),
        optimization_number = op
      )
    }
    
    map_correct_rot = Racmacs::rotateMap(
      map_rot,
      -(165-getMapAngle(map_rot, "Wu1", "XBB.1.5", op = op)),
      optimization_number = op
    )
    
    map$optimizations[[op]] = map_correct_rot$optimizations[[op]]
    
  }
  
  
  
  map
}
