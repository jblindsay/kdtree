import unittest
import random, times, strformat
import kdtree

test "kdtree tests":
  randomize()

  let numPoints = 10_000
  var
      points = newSeqOfCap[array[2, float]](numPoints)
      values = newSeqOfCap[int](numPoints)
      x: float
      y: float

  for a in 0..<numPoints:
      x = rand(100_000) / 1000
      y = rand(100_000) / 1000
      points.add([x, y])
      values.add(a)

  echo fmt"Building tree of {numPoints} random points..."
  var t0 = epochTime()
  var tree = newKdTree[int](points, values)
  var elapsed = epochTime() - t0
  echo fmt"Tree building elapsed time={elapsed}s"
  # echo fmt"Tree len={tree.len}"

  # perform nearestNeighour searches
  let numSearches = 10_000
  t0 = epochTime()
  for a in 0..<numSearches:
      x = rand(100_000) / 1000
      y = rand(100_000) / 1000
      discard tree.nearestNieghbour([x, y])
  
  elapsed = epochTime() - t0
  echo fmt"{numSearches} NN elapsed time={elapsed}s"

  t0 = epochTime()
  let k = 10
  for a in 0..<numSearches:
      x = rand(100_000) / 1000
      y = rand(100_000) / 1000
      discard tree.nearestNieghbours([x, y], k)

  elapsed = epochTime() - t0
  echo fmt"{numSearches} kNN (k={k}) elapsed time={elapsed}s"

  # import strabo/vector/shapefile
  # import strabo/primitives
  # let 
  #     testPointsFile = "/Users/johnlindsay/Documents/data/DigAg/yield/Farms/Roney 2015 Wheat/test_points.shp"
  #     hyperrectPointsFile = "/Users/johnlindsay/Documents/data/DigAg/yield/Farms/Roney 2015 Wheat/hyperrect.shp"
  #     inputFile = "/Users/johnlindsay/Documents/data/DigAg/yield/Farms/Roney 2015 Wheat/points_utm.shp"

  # echo "Reading points..."

  # var
  #     sf = openShapefile(inputFile)
  #     testPoints = openShapefile(testPointsFile)
  #     hyperrectPoints = openShapefile(hyperrectPointsFile)
  #     points = newSeqOfCap[array[K, float]](sf.numRecords)
  #     values = newSeqOfCap[int](sf.numRecords)
  #     point: Point

  # for p in 0..<sf.numRecords:
  #     if sf[p].shapeType != stNull:
  #         values.add(p+1)
  #         point = sf[p].points[0]
  #         points.add([point.x, point.y])

  # echo "Building tree..."
  # var tree = newKdTree[int](points, values)

  # echo fmt"tree height: {height(tree)}"
  # echo fmt"tree len: {len(tree)}"
  # for p in 0..<testPoints.numRecords:
  #     if testPoints[p].shapeType != stNull:
  #         point = testPoints[p].points[0]
  #         # tree.add([point.x, point.y], len(tree))
  #         echo fmt"fid={p+1}, x={point.x}, y={point.y}"

  #         # let (pt, values, dist) = tree.nearestNieghbour([point.x, point.y], squaredDist=false)
  #         # echo "Nearest Neighbour:"
  #         # echo fmt"point={pt}, value={values}, dist={dist}"

  #         # var ret = tree.nearestNieghbours([point.x, point.y], 4, squaredDist=false)
  #         # echo "Nearest 4 Neighbours:"
  #         # for (pt, value, dist) in ret:
  #         #     echo fmt"point={pt}, value={value}, dist={dist}"

  #         var ret2 = tree.withinRadius([point.x, point.y], radius=5.0, squaredDist=false, sortResults=true)
  #         echo "Within 5m Range:"
  #         for (pt, value, dist) in ret2:
  #             echo fmt"point={pt}, value={value}, dist={dist}"
          

  # var 
  #     hyperRects = newSeq[HyperRectangle]()
  #     min: array[K, float]
  #     max: array[K, float]
  #     hyperRect = 1
  # for p in 0..<hyperrectPoints.numRecords:
  #     if hyperrectPoints[p].shapeType != stNull:
  #         point = hyperrectPoints[p].points[0]
  #         tree.add([point.x, point.y], len(tree))
  #         if p mod 2 == 0:
  #             min[0] = point.x
  #             min[1] = point.y
  #         else:
  #             max[0] = point.x
  #             max[1] = point.y
  #             hyperRects.add(newHyperRectangle(min, max))
  #             var ret = tree.withinRange(hyperRects[hyperRect-1])
  #             echo fmt"Within HyperRectangle {hyperRect} Range:"
  #             for (pt, value) in ret:
  #                 echo fmt"point={pt}, value={value}"

  #             hyperRect += 1

  # echo fmt"Balance: {tree.isBalanced()}"
  # tree.rebalance()
  # echo fmt"After re-balance: {tree.isBalanced()}"

