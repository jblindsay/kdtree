# kdtree

**Contents**

1. [Description](#1-description)
2. [Documentation](#2-documentation)
3. [Usage](#2-usage)

## 1 Description

**kdtree** is a pure Nim [k-d tree](https://en.wikipedia.org/wiki/K-d_tree) implementation. k-d trees are data structures for performing efficient spatial query operations on point data sets. This implementation is very flexible, allowing for nearest-neighbour (single and multiple), within-radius (circular search areas), and range (rectangular search areas) spatial queries.

## 2 Documentation

Documentation for `kdtree` can be found [here](https://jblindsay.github.io/kdtree/kdtree.html).

## 3 Usage

```nim
import random, strformat
import kdtree
 
# Generate 100,000 random points
let numPoints = 100_000
var
  points = newSeqOfCap[array[2, float]](numPoints)
  values = newSeqOfCap[int](numPoints) 
  x: float
  y: float
  r = initRand(34)
 
for a in 0..<numPoints:
  x = r.rand(100.0)
  y = r.rand(100.0)
  points.add([x, y])
  values.add(a)
 
echo fmt"Building tree of {numPoints} random points..."
var tree = newKdTree[int](points, values)
# Notice that our 'values' are of int type here; the data associated with points can be of any generic data type.

# The preferred method of tree construction is bulk loading of point data using 'newKdTree'. However, you
# may also add individual points.
x = r.rand(100.0)
y = r.rand(100.0)
let value = numPoints
tree.add([x, y], value)

# However, adding many individual points can result in an unbalanced tree, which can result in inefficient queries. 
# You may check the tree balance and re-balance the tree if necessary.
let balance = tree.isBalanced() # The larger the value magnitude, the more unbalanced the tree is. The sign indicates 
                                # the direction of skew, with negative values indicating a left-skewed tree and positive 
                                # values indicated a right-skewed tree.

if abs(balance) > 1:
    tree.rebalance()


############################ 
# Spatial query operations #
############################ 

# Perform nearestNeighour searches
let numSearches = 10_000
for a in 0..<numSearches:
  x = r.rand(100.0)
  y = r.rand(100.0)
  let (pt, values, dist) = tree.nearestNeighbour([x, y])
  echo fmt"point={pt}, value={value}, dist={dist}"
 
# Perform nearestNeighours searches
let n = 10
for a in 0..<numSearches:
  x = r.rand(100.0)
  y = r.rand(100.0)
  let ret = tree.nearestNeighbours([x, y], n)
  for (pt, value, dist) in ret:
    echo fmt"point={pt}, value={value}, dist={dist}"

# Perform a withinRadius search
x = 50.0
y = 50.0
var ret2 = tree.withinRadius([x, y], radius=5.0, sortResults=true)
for (pt, value, dist) in ret2:
  echo fmt"point={pt}, value={value}, dist={dist}"
 
# Perform a withinRange search
var 
  min: array[2, float] = [0.0, 0.0]
  max: array[2, float] = [10.0, 10.0]
  hyperRect = newHyperRectangle(min, max)

var ret = tree.withinRange(hyperRect)
for (pt, value) in ret:
  echo fmt"point={pt}, value={value}"
