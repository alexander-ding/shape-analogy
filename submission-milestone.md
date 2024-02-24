## Mesh (milestone submission)

Please fill this out and submit your work to Gradescope by the milestone deadline.

### Mesh Validator

I'm using the halfedge representation, and I'm storing a linked list of all halfedges, vertices, faces, and edges.

The asserts are:

1. Every halfedge, vertex, edge, and face stored in the linked list should be a valid pointer.
1. Every halfedge, vertex, edge, and face backpoints to the correct position in the linked list.
1. Every halfedge should have a twin whose twin is the halfedge.
1. Every halfedge's twin should be different than the halfedge itself.
1. Every halfedge should link back to itself via exactly three next pointers (forming a triangle).
1. Every vertex's halfedge should point back to the vertex.
1. Every edge's halfedge should point to the edge.
1. Every edge's halfedge's twin should point to the edge.
1. Every face's halfedge should point back to the face.
1. Every face's halfedge's next halfedge should point back to the face.
1. Every face's halfedge's next next halfedge should point back to the face.
1. Every single pointer points to an address stored in the corresponding linked list (no untracked elements).

### Collaboration/References

N/A

### Known Bugs

N/A
