# Project 5: Motion Planning

This is the fifth project for CS 290: Computational Geometry at Duke University, fall 2024.

## Change Log
- no changes yet

## Table of Contents
- [Project 5: Motion Planning](#project-5-motion-planning)
  - [Change Log](#change-log)
  - [Table of Contents](#table-of-contents)
  - [Important Dates](#important-dates)
  - [Goals](#goals)
  - [System Requirements](#system-requirements)
  - [Organization](#organization)
    - [Helper Files](#helper-files)
  - [Figures](#figures)
  - [Submitting to Gradescope](#submitting-to-gradescope)
  - [Task 1: Vertical Ray Shooting (again) (10%)](#task-1-vertical-ray-shooting-again-10)
  - [Task 2: Vertical Decomposition Edges (70%)](#task-2-vertical-decomposition-edges-70)
  - [Task 4: Analysis (20%)](#task-4-analysis-20)


## Important Dates
- Tuesday 11/26: Project released
- Sunday 12/1 @ 11:59pm: Progress report **due** (make a private post on Ed with subject "Check-In")
- Friday 12/6 @ 11:59pm: Project **due**

## Goals
For this project you will complete an implementation of motion planning for a ``point'' robot among polygonal obstacles, contained in a rectangular region. In other words, the workspace is a rectangle with polygonal holes. See sections 6.1 and 13.2 of the textbook for more information, however the project can be completed without any of the implementation details therein as we provided much of them for you.

## System Requirements

This project leverages more packages than in the past. If you are using the Anaconda distribution then these should already be installed. If not, you can install them with `pip`:

`pip install matplotlib networkx`

Depending on how your Python installation is set up on your system, you may need to ensure that `pip` has been installed and is on your PATH.

<!-- ## Object-Oriented Python
The project heavily uses (basic) object-oriented principles (OOP) to organize the flow of logic. For one of many excellent resources for learning about Python's OOP design, we recommend the [official Python tutorial](https://docs.python.org/3/tutorial/classes.html), but there are many others. -->

## Organization
This repository includes all necessary skeleton files, as well as some example snippets that may be useful for testing your code and/or implementing the following tasks. In particular, we have provided some code for basic uses of `matplotlib` that you may use as starting points for your own visualization code. All parts use the `primitives.py` file and the classes defined therein.

### Helper Files

- `primitives.py`: Contains `Point`, `Segment`, `Line` classes, as seen in previous projects. These objects support arbitrary-precision arithmetic **if they have rational coordinates**, which all of our tests and provided code use.  All classes provide a `.draw()` method for visualization purposes.
  
## Figures

![](./figs/workspace.png)

Figure 1: An example workspace with two polygonal obstacles, which are represented by two inner cycles of the polygon.

![](./figs/decomp.png)

Figure 2: The vertical decomposition of the workspace, where the red segments erect downwards from obstacle vertices and blue segments erect upwards from obstacle vertices. Note that there are no vertical segments inside the inner cycles.

![](./figs/roadmap.png)

Figure 3: The roadmap, where we have a vertex per trapezoid and each vertical edge of the trapezoids.

![](./figs/path.png)

Figure 4: The shortest path, by number of edges, from $s$ to $t$ in the roadmap graph. The first and last edges are to/from the centers of the trapezoids that contain $s,t$.

## Submitting to Gradescope
For this Project we will have two "assignments" on Gradescope: one for all project code, and one to upload your PDF (see Task 4 - Analysis). You should upload **one** file, `motion_planning.py`, to the Code portion.

## Task 1: Vertical Ray Shooting (again) (10%)

For this task, complete the method `naive_vertical_shoot()` in `motion_planning.py`. This method has *minor differences* from the vertical ray shooting method from Project 4. Specifically, one should return four items, the first two being the visible segment and the point on it upwards from the given query point, and the last two being the visible segment and the point on it downards from the given query point. **Furthermore, we do not consider a segment containing a query point to be visible**, otherwise this method would trivially always return the segment that a query vertex belongs to.

## Task 2: Vertical Decomposition Edges (70%)

For this task, complete the method `all_vd_queries()` in `motion_planning.py`. Given the outer cycle `outer` and list of inner cycles `inners`, return a representation of the vertical decomposition edges. Namely, for each such edge, represent it as a 3-tuple: `(q, visible_seg, visible_point)` where `q` is a vertex of an inner cycle, `visible_seg` is the visible segment above or below `q`, and `visible_point` is the point on that segment visible from `q`. **However, we do not want any edges that lie inside any inner cycle**, so do not include such corresponding tuples in your output. To determine if, say, the segment emanating upwards from `q` to its upwards-visible segment lies on the same inner cycle as `q` (and thus the segment lies in that inner cycle), it suffices to look at how the two edges of the inner cycle appear around `q` in relation to the vertical segment. This check is expected to run in only constant time.

For a visualization, please see the red and blue segments of the vertical decomposition in Figure 2 of the [Figures section](#figures).

## Task 4: Analysis (20%)

For the analysis section, we simply ask for figures similar to those seen in the [Figures section](#figures), but for the example input described at the bottom of `motion_planning.py`. With a complete implementation, running the script will create each of these four figures, one-by-one.

1) **Figure/Caption**: The example workspace as defined by `outer` and `inners`.
2) **Figure/Caption**: The vertical decomposition.
3) **Figure/Caption**: The roadmap.
4) **Figure/Caption**: The path from $s = (0.5,6)$ and $t = (7,8.5)$.)