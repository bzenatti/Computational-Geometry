# Computational Geometry Exercises

This repository contains exercises and projects focused on computational geometry. It is developed using an adapted version of the CGAL (Computational Geometry Algorithms Library) by Ricardo Dutra da Silva at the Federal University of Technology – Paraná, Brazil.

## Overview

The repository offers hands-on experience with implementing algorithms and solving problems in computational geometry. It includes a variety of exercises that cover key concepts, as well as two main projects that involve triangulation and mesh manipulation.

## Project: Triangulating and Approximating an Image

- **Image Conversion:**  
  The project reads an input image, converts it to grayscale (in PGM format), and resizes it to 512×512 pixels (if needed) using OpenCV.

- **Point Selection:**  
  A specified number of sample points are selected from the image (with an option to include border points). These points serve as the vertices for the triangulation.

- **Delaunay Triangulation:**  
  - Constructs a super triangle that encloses all the sample points.
  - Incrementally inserts points into the mesh, either inside a triangle or along an edge.
  - Subdivides triangles and legalizes edges (using in-circle tests and edge flips) to maintain the Delaunay condition.
  - Removes the super triangle and adds border edges to finalize the mesh.

- **Mesh Rasterization:**  
  The triangulated mesh is rasterized by interpolating pixel intensities using barycentric coordinates, resulting in an approximated version of the original image.
