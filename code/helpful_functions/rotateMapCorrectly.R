rotateMapCorrectly = function(map){
  
  getMapAngle = function(map, ag1, ag2){
    coords = Racmacs::agCoords(map)[c(ag1, ag2),]
    x = atan(diff(coords[,2]) / diff(coords[,1])) * 180 / pi
    if (coords[2,1]> coords[1,1]) x = 180 + x
    x
  }
  
  map_rot = Racmacs::rotateMap(
    map,
    -(180-getMapAngle(map, "Wu-1", "XBB.1.5"))
  )
  
  if (agCoords(map_rot)["Wu-1",2] < agCoords(map_rot)["JN.1",2]){
    map_ref = Racmacs::reflectMap(map)
    map_rot = Racmacs::rotateMap(
      map_ref,
      -(180-getMapAngle(map_ref, "Wu-1", "XBB.1.5"))
    )
  }
  
  map_correct_rot = Racmacs::rotateMap(
    map_rot,
    -(165-getMapAngle(map_rot, "Wu-1", "XBB.1.5"))
  )
  
  map_correct_rot
}