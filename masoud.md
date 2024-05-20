## New classes and files
Update: [20.05.2024]
I have used std::vector instead of pointers to arrays, and std::array instead of fixed arrays. 
###  namespace utils
1. enums for program options: program_modes, accept_modes, polarization_modes, program_options
2. Simple random generators in c++ style
3. Constants from the old const.h
4. absolute_error and relative_error templates
5. read_cmd: reads the command arguments
6. show_progress: shows a simple progress bar for long operations
7. linspace: generates a range for P_t, y_p and so on from Andrea
8. four_vec simple four vector std::array, and r2_tensor simple rank 2 Minkowski tensor. All the tensors are assumed to have lower indices. Some required functions are added.
## namespace  utils::geometry
Here I defined a new class four_vector which encapsulates a Minkowski four vector with keeping the index structure in check. It contains the array utils::four_vec as internal data. Some important points:
1. I optimized the vector addition without going too far, such as using [expression templates](https://en.wikipedia.org/wiki/Expression_templates)  which is an overkill. 
[![](https://images.unsplash.com/photo-1493612276216-ee3925520721?q=80&w=1000&auto=format&fit=crop&ixlib=rb-4.0.3&ixid=M3wxMjA3fDB8MHxzZWFyY2h8Mnx8cmFuZG9tfGVufDB8fDB8fHww)](https://images.unsplash.com/photo-1493612276216-ee3925520721?q=80&w=1000&auto=format&fit=crop&ixlib=rb-4.0.3&ixid=M3wxMjA3fDB8MHxzZWFyY2h8Mnx8cmFuZG9tfGVufDB8fDB8fHww)

![utils benchmark][img/benchutils.png]

