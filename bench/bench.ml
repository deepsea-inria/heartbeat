open XBase
open Params

let system = XSys.command_must_succeed_or_virtual

(*****************************************************************************)
(** Parameters *)

let arg_virtual_run = XCmd.mem_flag "virtual_run"
let arg_virtual_build = XCmd.mem_flag "virtual_build"
let arg_virtual_get = XCmd.mem_flag "virtual_get"
let arg_force_get = XCmd.mem_flag "force_get"
let arg_nb_runs = XCmd.parse_or_default_int "runs" 1
let arg_mode = Mk_runs.mode_from_command_line "mode"
let arg_skips = XCmd.parse_or_default_list_string "skip" []
let arg_onlys = XCmd.parse_or_default_list_string "only" []
let arg_show_utilization = false
let arg_benchmarks = XCmd.parse_or_default_list_string "benchmark" ["all"]
let arg_proc =
  let hostname = Unix.gethostname () in
  let default =
    if hostname = "teraram" then
      [ 40; ]
    else if hostname = "cadmium" then
      [ 48; ]
    else if hostname = "hiphi.aladdin.cs.cmu.edu" then
      [ 64; ]
    else if hostname = "aware.aladdin.cs.cmu.edu" then
      [ 72; ]
    else if hostname = "beast" then
      [ 8; ]
    else
      [ 1; ]
  in
  let default =
    if List.exists (fun p -> p = 1) default then
      default
    else
      1 :: default
  in
  XCmd.parse_or_default_list_int "proc" default
let arg_print_err = XCmd.parse_or_default_bool "print_error" false
let arg_scheduler = XCmd.parse_or_default_string "scheduler" ""
let arg_nb_make_cores =
  let proc = List.fold_left max 1 arg_proc in
  XCmd.parse_or_default_int "nb_make_cores" proc
let arg_path_to_data = XCmd.parse_or_default_string "path_to_data" "_data"
let arg_path_to_results = XCmd.parse_or_default_string "path_to_results" "_results"
                
let run_modes =
  Mk_runs.([
    Mode arg_mode;
    Virtual arg_virtual_run;
    Runs arg_nb_runs; ])

let multi_proc = List.filter (fun p -> p <> 1) arg_proc

(*****************************************************************************)
(** Steps *)

let select get make run check plot =
   let arg_skips =
      if List.mem "run" arg_skips && not (List.mem "make" arg_skips)
         then "make"::arg_skips
         else arg_skips
      in
   let run' () = (
       system (Printf.sprintf "mkdir -p %s" arg_path_to_results) false;
       run())
   in
   Pbench.execute_from_only_skip arg_onlys arg_skips [
      "get", get;
      "make", make;
      "run", run';
      "check", check;
      "plot", plot;
      ]

let nothing () = ()

(*****************************************************************************)
(** Files and binaries *)

let build path bs is_virtual =
   system (sprintf "make -C %s -j%d %s" path arg_nb_make_cores (String.concat " " bs)) is_virtual

let file_results exp_name =
  Printf.sprintf "%s/results_%s.txt" arg_path_to_results exp_name

let file_tables_src exp_name =
  Printf.sprintf "%s/tables_%s.tex" arg_path_to_results exp_name

let file_tables exp_name =
  Printf.sprintf "%s/tables_%s.pdf" arg_path_to_results exp_name

let file_plots exp_name =
  Printf.sprintf "%s/plots_%s.pdf" arg_path_to_results exp_name

(** Evaluation functions *)

let eval_exectime = fun env all_results results ->
  Results.get_mean_of "exectime" results

let eval_exectime_stddev = fun env all_results results ->
  Results.get_stddev_of "exectime" results

let string_of_millions ?(munit=false) v =
   let x = v /. 1000000. in
   let f = 
     if x >= 10. then sprintf "%.0f" x
     else if x >= 1. then sprintf "%.1f" x
     else if x >= 0.1 then sprintf "%.2f" x
     else sprintf "%.3f" x in
   f ^ (if munit then "m" else "")
                        
let formatter_settings = Env.(
    ["prog", Format_custom (fun s -> "")]
  @ ["algorithm", Format_custom (fun s -> s)]
  @ ["n", Format_custom (fun s -> sprintf "Input: %s million 32-bit ints" (string_of_millions (float_of_string s)))]
  @ ["proc", Format_custom (fun s -> sprintf "#CPUs %s" s)]
  @ ["promotion_threshold", Format_custom (fun s -> sprintf "F=%s" s)]
  @ ["threshold", Format_custom (fun s -> sprintf "K=%s" s)]
  @ ["block_size", Format_custom (fun s -> sprintf "B=%s" s)]      
  @ ["operation", Format_custom (fun s -> s)])

let default_formatter =
  Env.format formatter_settings
    
let string_of_percentage_value v =
  let x = 100. *. v in
  (* let sx = if abs_float x < 10. then (sprintf "%.1f" x) else (sprintf "%.0f" x)  in *)
  let sx = sprintf "%.1f" x in
  sx
    
let string_of_percentage ?(show_plus=true) v =
   match classify_float v with
   | FP_subnormal | FP_zero | FP_normal ->
       sprintf "%s%s%s"  (if v > 0. && show_plus then "+" else "") (string_of_percentage_value v) "\\%"
   | FP_infinite -> "$+\\infty$"
   | FP_nan -> "na"

let string_of_percentage_change ?(show_plus=true) vold vnew =
  string_of_percentage ~show_plus:show_plus (vnew /. vold -. 1.0)

let ipfs_get hash outfile is_virtual =
  system (sprintf "wget -O %s https://ipfs.io/ipfs/%s" outfile hash) is_virtual

let ipfs_get_if_needed hash outfile force_get is_virtual =
  if force_get || not (Sys.file_exists outfile) then
    ipfs_get hash outfile is_virtual
  else
    ()

let ipfs_get_files table force_get is_virtual =
  List.iter (fun (h, p) -> ipfs_get_if_needed h p force_get is_virtual) table

(*****************************************************************************)
(** Sequence-library benchmark *)

module ExpSequenceLibrary = struct

let name = "sequence"

let prog_heartbeat = name^".heartbeat"

let prog_cilk = name^".cilk"

let mk_pbbs_algorithm = mk_prog prog_cilk & mk string "algorithm" "pbbs"
                           
let mk_algorithms = (
(*     (mk_prog prog_heartbeat & mk string "algorithm" "sequential")
  ++ *) mk_pbbs_algorithm
(*  ++ (  (mk_prog prog_heartbeat & mk string "algorithm" "heartbeat_sequential")
      & mk_heartbeat_settings) *)
  ++ (  (mk_prog prog_heartbeat & mk string "algorithm" "heartbeat"))
)

let mk_operations = mk_list string "operation" [ "reduce"; "max_index"; "scan"; "pack"; "filter"; ] 

let mk_input_sizes = mk_list int "n" [ 400000000 ]

let mk_procs = mk_list int "proc" [ 1; (* 10; 20; 30; *) 40; ]

let get() = ()

let make() =
  build "." [prog_heartbeat; prog_cilk] arg_virtual_build

let run() =
  Mk_runs.(call (run_modes @ [
    Output (file_results name);
    Timeout 400;
    Args (mk_algorithms & mk_procs & mk_operations & mk_input_sizes)
  ]))

let check = nothing  (* do something here *)

let plot() =
     Mk_bar_plot.(call ([
      Bar_plot_opt Bar_plot.([
                              (* Chart_opt Chart.([Dimensions (12.,8.) ]);*)
                              X_titles_dir Vertical;
         Y_axis [ Axis.Lower (Some 0.); Axis.Upper (Some 2.5);
                  Axis.Is_log false ] 
         ]);
      Formatter default_formatter;
      Charts (mk_procs & mk_input_sizes);
      Series mk_algorithms;
      X mk_operations;
      Y_label "Time (s)";
      Y eval_exectime;
      Y_whiskers eval_exectime_stddev;
      Output (file_plots name);
      Results (Results.from_file (file_results name));
      ]))

let all () = select get make run check plot

end

(*****************************************************************************)
(** Comparison benchmark *)

module ExpCompare = struct

let name = "compare"

let all_benchmarks =
  match arg_benchmarks with
  | ["all"] -> [
    "convexhull"; "samplesort"; "radixsort"; "nearestneighbors";
    "suffixarray"; "removeduplicates"; (*"mis";*) "mst"; "matching"; "spanning";
    "delaunay"; (*"bfs";*) (*"refine"; *) "raycast"; (*"pbfs";*)
    ]
  | _ -> arg_benchmarks
    
let heartbeat_prog_of n = n ^ ".heartbeat"
let cilk_prog_of n = n ^ ".cilk"
let cilk_elision_prog_of n = n ^ ".cilk_elision"
                             
let heartbeat_progs = List.map heartbeat_prog_of all_benchmarks
let cilk_progs = List.map cilk_prog_of all_benchmarks
let cilk_elision_progs = List.map cilk_elision_prog_of all_benchmarks
let all_progs = List.concat [heartbeat_progs; cilk_progs; cilk_elision_progs]

let path_to_infile n = arg_path_to_data ^ "/" ^ n

let mk_infiles ty descr = fun e ->
  let f (p, t, n) =
    let e0 = 
      Env.add Env.empty "infile" (string p)
    in
    Env.add e0 ty t
  in
  List.map f descr

(*****************)
(* Convex hull *)
    
let input_descriptor_hull = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "array_point2d_in_circle_large.bin", string "2d", "in circle";
  "array_point2d_kuzmin_large.bin", string "2d", "kuzmin";
  "array_point2d_on_circle_medium.bin", string  "2d", "on circle";
]

let mk_hull_infiles = mk_infiles "type" input_descriptor_hull
                                 
let mk_heartbeat_prog n =
  let a =
    (mk string "prog" (heartbeat_prog_of n))
      & (mk string "algorithm" "heartbeat")
  in
  if arg_scheduler = "" then
    a
  else
    a & (mk string "scheduler" arg_scheduler)

let mk_pbbs_lib =
  mk string "algorithm" "pbbs"
    
let mk_pbbs_prog n =
    (mk string "prog" (cilk_prog_of n))
  & mk_pbbs_lib

let mk_single_proc = mk int "proc" 1

let mk_multi_proc = mk_list int "proc" multi_proc

let mk_pbbs_elision_prog n =
    (mk string "prog" (cilk_elision_prog_of n))
  & mk_pbbs_lib
    
type input_descriptor =
    string * Env.value * string (* file name, type, pretty name *)
    
type benchmark_descriptor = {
  bd_name : string;
  bd_infiles : Params.t;
  bd_input_descr : input_descriptor list;
}

(*****************)
(* Sample sort *)
  
let input_descriptor_samplesort = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "array_double_random_large.bin", string "double", "random";
  "array_double_exponential_large.bin", string "double", "exponential";
  "array_double_almost_sorted_10000_large.bin", string "double", "almost sorted";
]
    
let mk_samplesort_infiles = mk_infiles "type" input_descriptor_samplesort

(*****************)
(* Radix sort *)

let input_descriptor_radixsort = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "array_int_random_large.bin", string "int", "random";    
  "array_int_exponential_large.bin", string "int", "exponential";
  "array_pair_int_int_random_256_large.bin", string "pair_int_int", "random pair" (*"random int pair 256" *);
(*  "array_pair_int_int_random_100000000_large.bin", string "pair_int_int", "random int pair 10m";*)
]

let mk_radixsort_infiles = mk_infiles "type" input_descriptor_radixsort

(*****************)
(* BFS *)
      (*
let input_descriptor_bfs = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "cube_large.bin", int 0, "cube";
  "rmat24_large.bin", int 0, "rMat24";
(*  "rmat27_large.bin", int 0, "rMat27";*)
]

let mk_bfs_infiles = mk_infiles "source" input_descriptor_bfs
    *)
(*****************)
(* PBFS *)
(*
let input_descriptor_pbfs = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "cube_large.bin", int 0, "cube";
  "rmat24_large.bin", int 0, "rMat24";
(*  "rmat27_large.bin", int 0, "rMat27";*)
]

let mk_pbfs_infiles = mk_infiles "source" input_descriptor_pbfs
*)
(*****************)
(* MIS *)

let input_descriptor_mis = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "cube_large.bin", int 0, "cube";
  "rmat24_large.bin", int 0, "rMat24";
(*  "rmat27_large.bin", int 0, "rMat27";*)
]

let mk_mis_infiles = mk_infiles "source" input_descriptor_mis

(*****************)
(* MST *)

let input_descriptor_mst = input_descriptor_mis

let mk_mst_infiles = mk_infiles "source" input_descriptor_mst

(*****************)
(* Matching *)
(*
let input_descriptor_matching = input_descriptor_mis

let mk_matching_infiles = mk_infiles "source" input_descriptor_matching
*)
(*****************)
(* Spanning *)

let input_descriptor_spanning = input_descriptor_mis

let mk_spanning_infiles = mk_infiles "source" input_descriptor_spanning

(*****************)
(* Suffix array *)

let input_descriptor_suffixarray = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "chr22.dna.bin", string "string", "dna";
  "etext99.bin", string "string", "etext";
  "wikisamp.xml.bin", string "string", "wikisamp";
]
      
let mk_suffixarray_infiles = mk_infiles "type" input_descriptor_suffixarray

(*****************)
(* Nearest neighbors *)

let input_descriptor_nearestneighbors = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "array_point2d_kuzmin_medium.bin", string "array_point2d", "kuzmin";
(*  "array_point3d_on_sphere_medium.bin", string "array_point3d", "on sphere";*)
  "array_point3d_plummer_medium.bin", string "array_point3d", "plummer"; 
  (*  "array_point2d_in_square_medium.bin", string "array_point2d", "in square";*)
(*  "array_point3d_in_cube_medium.bin", string "array_point3d", "in cube"; *)
]

let mk_nearestneighbors_infiles = mk_infiles "type" input_descriptor_nearestneighbors

(*****************)
(* Delaunay *)

let input_descriptor_delaunay = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "array_point2d_in_square_delaunay_large.bin", string "array_point2d", "in square";
  "array_point2d_kuzmin_delaunay_large.bin", string "array_point2d", "kuzmin";
]

let mk_delaunay_infiles = mk_infiles "type" input_descriptor_delaunay

(*****************)
(* Refine *)

let input_descriptor_refine = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "triangles_point2d_delaunay_in_square_refine_large.bin", string "triangles_point2d", "in square";
  "triangles_point2d_delaunay_kuzmin_refine_large.bin", string "triangles_point2d", "kuzmin";
]

let mk_refine_infiles = mk_infiles "type" input_descriptor_refine    

(*****************)
(* Raycast *)

let input_descriptor_raycast = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "happy_ray_cast_dataset.bin", string "raycast", "happy";
  "xyzrgb_manuscript_ray_cast_dataset.bin", string "raycast", "xyzrgb";
]

let mk_raycast_infiles = mk_infiles "type" input_descriptor_raycast    

(*****************)
(* Remove duplicates *)

let input_descriptor_removeduplicates = List.map (fun (p, t, n) -> (path_to_infile p, t, n)) [
  "array_int_random_large.bin", string "array_int", "random";
  "array_int_random_bounded_100000_large.bin", string "array_int", "bounded random";
  "array_int_exponential_large.bin", string "array_int", "exponential";
  "array_string_trigrams_large.bin", string "array_string", "string trigrams";
]

let mk_removeduplicates_infiles = mk_infiles "type" input_descriptor_removeduplicates

(*****************)
(* All benchmarks *)

let benchmarks' : benchmark_descriptor list = [
  { bd_name = "radixsort";
    bd_infiles = mk_radixsort_infiles;
    bd_input_descr = input_descriptor_radixsort;
  };
  { bd_name = "samplesort";
    bd_infiles = mk_samplesort_infiles;
    bd_input_descr = input_descriptor_samplesort;
  };
  { bd_name = "suffixarray";
    bd_infiles = mk_suffixarray_infiles;
    bd_input_descr = input_descriptor_suffixarray;
  };
  { bd_name = "removeduplicates";
    bd_infiles = mk_removeduplicates_infiles;
    bd_input_descr = input_descriptor_removeduplicates;
  };
  { bd_name = "convexhull";
    bd_infiles = mk_hull_infiles;
    bd_input_descr = input_descriptor_hull;
  };
  { bd_name = "nearestneighbors";
    bd_infiles = mk_nearestneighbors_infiles;
    bd_input_descr = input_descriptor_nearestneighbors;
  };
  { bd_name = "delaunay";
    bd_infiles = mk_delaunay_infiles;
    bd_input_descr = input_descriptor_delaunay;
  };
  { bd_name = "refine";
    bd_infiles = mk_refine_infiles;
    bd_input_descr = input_descriptor_refine;
  };
  { bd_name = "raycast";
    bd_infiles = mk_raycast_infiles;
    bd_input_descr = input_descriptor_raycast;
  }; (*
  { bd_name = "bfs";
    bd_infiles = mk_bfs_infiles;
    bd_input_descr = input_descriptor_bfs;
  }; 
  { bd_name = "pbfs";
    bd_infiles = mk_pbfs_infiles;
    bd_input_descr = input_descriptor_pbfs;
  }; *)
  { bd_name = "mis";
    bd_infiles = mk_mis_infiles;
    bd_input_descr = input_descriptor_mis;
  }; 
  { bd_name = "mst";
    bd_infiles = mk_mst_infiles;
    bd_input_descr = input_descriptor_mst;
  }; (*
  { bd_name = "matching";
    bd_infiles = mk_matching_infiles;
    bd_input_descr = input_descriptor_matching;
  }; *)
  { bd_name = "spanning";
    bd_infiles = mk_spanning_infiles;
    bd_input_descr = input_descriptor_spanning;
  }; 
  
]

let infiles_by_hash = [
  "QmXZjB1y8uFZ5RjwsiA9JvjyoCNBHwKAvKFtGY7rb7tA5V", "3Dgrid_J_10000000.bin";
  "QmZFTC6Zbi9qyJLAyprPc6z8GghvkDcNis9Sta2QrPjX1j", "angel_ray_cast_dataset.bin";
  "QmXv2RnFr1H5S4ip3LQoZxQffL89LUjmrc93ehui21bkUw", "array_double_almost_sorted_10000_large.bin";
  "QmYLMFaXDKvz7kS5uUa1CHMKi1kYCpiUhNzvCQKakD9ei1", "array_double_almost_sorted_1000_small.bin";
  "QmbuUCyLXNrF15uCVaerHxKNXXDfYknSt76ff1SZyw9N9Y", "array_double_almost_sorted_3162_medium.bin";
  "Qmb3RJMbwQir2mVKYYSA3wzmD7ZeXTUFUXt7VDo5xkErKe", "array_double_exponential_large.bin";
  "QmY5Rd9ovjc5NC5aZEQm6XMTsFPXuK5mcETQUtU3YT7vSv", "array_double_exponential_medium.bin";
  "QmayPLeexUXd5CDU1toFEKD3GEN2RarR8rM6SYqETtGW2E", "array_double_exponential_small.bin";
  "QmQAxEYwvdDeTtU6ReoGUSpJCgJ477WvCHs7YqayFKMjcm", "array_double_random_large.bin";
  "QmabSxwMqL98tMZeaqBRxGUxyMVEpnq8ez8VCxQEp2geNb", "array_double_random_medium.bin";
  "QmQoMHANBS4zEZZEDMmRcD8dkdtLmtuuaqCcBtmyfAEa9x", "array_double_random_small.bin";
  "QmVDwFLa3USQMSVv3nfUBwjbsfo1SmATw9okjJT5AWrb4n", "array_int_exponential_large.bin";
  "QmNyPhMzvF54BojTuGqrQEaJb9XSmc6bpq2TrtNa3wpr7u", "array_int_exponential_medium.bin";
  "QmX5xninsAtQkJSSWJFWNF9yqwmDhmF9oCxFDKsia3jcrf", "array_int_exponential_small.bin";
  "QmcZ42Suo1AynTmmB8zum8VJ3FNhb59gMi4VyA6tYWe6gE", "array_int_random_bounded_100000_large.bin";
  "QmZY811Agj461by2Y1dqPuRaW8FeCDpZ62SwpxJub99KcE", "array_int_random_bounded_100000_medium.bin";
  "QmZyk1owLgLBQNELsgN7c14rXqJzAyj1jyB2QjdtvFdpok", "array_int_random_bounded_100000_small.bin";
  "QmUhvPaavdgMWGzFTpXRMaDwdRzNSrVSbdLHpMZsMnkrwG", "array_int_random_large.bin";
  "QmT8c1EE9GPurEf1P3gHoBJParyK4dWs1nQ1E8r2BpXcXm", "array_int_random_medium.bin";
  "QmPYB5cBdVLq4t84jU1MYcM5rZD3vwtqGU2uxmndUE7Jmv", "array_int_random_small.bin";
  "QmaMkDBQhTywM1t1QwLwPPDE336fBBzSB7QLY4UYsjdUpV", "array_pair_int_int_random_100000000_large.bin";
  "QmbxJ3Xny3yK1N11oKutAVELsnMYpKfAT31eQ7J8Emujiq", "array_pair_int_int_random_10000000_medium.bin";
  "QmbRrXoVjvuNznQ2Y8GbyZXs6ALRhX6y7ErHiGn1YPqAJ7", "array_pair_int_int_random_1000000_small.bin";
  "QmZt4u8McsZFdkpg9jZnvD6bcRitoPCCXARXuKX1YA9efV", "array_pair_int_int_random_256_large.bin";
  "QmSV3ofTHXJsHwkKAHabexZUdBPEXbDq7qTqRpDfhCZXuM", "array_pair_int_int_random_256_medium.bin";
  "QmXxWH1qwYDuHNTyaTM1b65bzFruK6tc7edGYmHhfRtao9", "array_pair_int_int_random_256_small.bin";
  "QmbBqaM5SAuuhx7YAY1PeyH62GbbiHWfPxYMmmiQdw22Pg", "array_pair_string_int_trigrams_large.bin";
  "QmUvgq6Vo6Mvo3vgXvxzHY8BZQvLmH1hBxEXzT3ub6PZaE", "array_pair_string_int_trigrams_medium.bin";
  "QmPuJfUaeBWiwL71kw91fgNNFJx1wWfvfgR875sugdbd3e", "array_pair_string_int_trigrams_small.bin";
  "QmRbzUS4eaVyBrNMtSi2xk5bgu8Ka9xTcti1Xi2aVwz3oK", "array_point2d_in_circle_large.bin";
  "QmP5g6Qwxuw34paLrvctLYR4GWvntraTaFwfqhbF2rfLtC", "array_point2d_in_circle_medium.bin";
  "QmfWVA6G9HhfLSivV9nPF2kZLCqsW8hiViLkvaK7idjP1K", "array_point2d_in_circle_small.bin";
  "QmaA4Jq23PEPL3dyv7F6soNQR15PDW8NhFdDZRzmuJapzJ", "array_point2d_in_square_delaunay_large.bin";
  "QmaA4Jq23PEPL3dyv7F6soNQR15PDW8NhFdDZRzmuJapzJ", "array_point2d_in_square_delaunay_medium.bin";
  "QmPYZumkzCVUeATptYKMNv5qEFgRpTPGYtZ35xkhc7BpJX", "array_point2d_in_square_delaunay_small.bin";
  "QmNex8BfEQweJPeposH1NS8KkJ4UuKg17FvNvvmDN34oNn", "array_point2d_in_square_large.bin";
  "QmaA4Jq23PEPL3dyv7F6soNQR15PDW8NhFdDZRzmuJapzJ", "array_point2d_in_square_medium.bin";
  "QmPYZumkzCVUeATptYKMNv5qEFgRpTPGYtZ35xkhc7BpJX", "array_point2d_in_square_small.bin";
  "QmSLAWyJ3kmFBmXSV7HvyRfiLgwkwfQm1JKKSYUMpQzpTi", "array_point2d_kuzmin_delaunay_large.bin";
  "QmSLAWyJ3kmFBmXSV7HvyRfiLgwkwfQm1JKKSYUMpQzpTi", "array_point2d_kuzmin_delaunay_medium.bin";
  "QmXYtFbsG6KRXmAZyuFiPaUKAQE7EQKH4y8gVMZg1KCACF", "array_point2d_kuzmin_delaunay_small.bin";
  "QmZ9UTmjSmPcEFxmH441fqx6g1LLej4rpVWFYtKbzyZLq9", "array_point2d_kuzmin_large.bin";
  "QmSLAWyJ3kmFBmXSV7HvyRfiLgwkwfQm1JKKSYUMpQzpTi", "array_point2d_kuzmin_medium.bin";
  "QmXYtFbsG6KRXmAZyuFiPaUKAQE7EQKH4y8gVMZg1KCACF", "array_point2d_kuzmin_small.bin";
  "Qma1oD5ojjSMgLdywLPzLJ5cBV8esUhcACWmN5SyUr9hFT", "array_point2d_on_circle_large.bin";
  "Qmcz9JbR9piozbTwugioAymyb7w7HjybYrAhDss3n6QkMv", "array_point2d_on_circle_medium.bin";
  "QmXJ3qDNgaFW6S2pFz8Td2ETVjUtDe8aaxDGYfeT2aQSsA", "array_point2d_on_circle_small.bin";
  "QmZrL33umj2k31CpRW1yG2YkeWfQue5cEN8YvXQqLJzTBe", "array_point3d_in_cube_large.bin";
  "QmVZRcB1CxpmTYq5ChCkWjRNYQWotyNvhCEZvmjmq5bHQS", "array_point3d_in_cube_medium.bin";
  "QmZDDGg4ucoo8znhiewJQSdJvdhMbRtEF3wA17Racq6Xtn", "array_point3d_in_cube_small.bin";
  "QmRDNcgGWVZyeauDSd867CTXo8V9ZuEp7K8DuBpmJY5duF", "array_point3d_on_sphere_large.bin";
  "QmT183mJMpsrVvNY9r4Pdyn7vTJErXyizj9BLh4TKgq3bj", "array_point3d_on_sphere_medium.bin";
  "QmZbp3ZWuP1DuXBsUV6uR1WhNiBdr2xUPoxKYnJ3rfVWQT", "array_point3d_on_sphere_small.bin";
  "QmX6W5b8F5kae38Muy9yxaEfsiqb6gHhQgiQa1VgVt9jKQ", "array_point3d_plummer_large.bin";
  "QmVZsRyDFRNECsxbrogSUdX5VrQwXCj7AMZPzizwuW1x4f", "array_point3d_plummer_medium.bin";
  "Qmesw1McsjPdkENNYomTUsTbMrVKpF3NQec5tFtUBk6bvB", "array_point3d_plummer_small.bin";
  "QmX4wdLHgpoGv92SXerXkgLqjZVF4gZYcnjtAbS51T4unP", "array_string_trigrams_large.bin";
  "QmPeRz42Axz5V1Jzyutet3YoDMTNdAXQWtyWGDcU6rRs6J", "array_string_trigrams_medium.bin";
  "QmacvDqCeAdQ4bZemHnXsN8pAPUVkK3JMS2SND6zihajiN", "array_string_trigrams_small.bin";
  "Qma4z37vrhKTiAXBUnaUeJS9cfrD6bJ256yHuGsfft1M5m", "chain_large.bin";
  "QmUGQqKtivcCYEQxZdEow2GJT3CgfkiijgcLUxdUKeJSYt", "chain_small.bin";
  "QmXVka2FKr6vHS5h8sr2L6wPjyD72EsQDbBDFghLn7sj9M", "chr22.dna.bin";
  "QmcACNqugJsNXrHK3wrsGedEV44tHEtdZXEn3zY1Eo9Qdj", "cube_large.bin";
  "QmfYtSDcXygyErVG5YfM4vXzLfAtr2k3BF9L7D1zbBFcgu", "cube_medium.bin";
  "QmPDhVjgQQG7FJExGPyRdSmRPhbcdT1TyVzvXgedBWnKjX", "cube_small.bin";
  "QmbpWbvT1zPXtwBCgxh7d6SdLV9YFeo9KTNNUqdSFyNfZR", "delaunay.bin";
  "QmcnrU8LhBg8P3uGBv5DzG1SWkQwihnMB7dDTHs6p419AE", "dragon_ray_cast_dataset.bin";
  "QmfX2ZG8TqqzDjYMUWXFKQujzXhcc4hzQZqmuTaHF9eT8t", "etext99.bin";
  "QmbEJ3QpiVQDWBwrnLCFBhkNQMcRJSKrs8TaguCnYUCR4r", "europe.bin";
  "QmaroxYoU5sJBVuhmZ51UttBfJpLuZ7yUDpXqNapfoUSkd", "grid_sq_large.bin";
  "QmRDaaCrsRji15z2f7PhJTJpvS78L7MEtSuAnhtk5aMxHm", "happy_ray_cast_dataset.bin";
  "QmRDaaCrsRji15z2f7PhJTJpvS78L7MEtSuAnhtk5aMxHm", "happy_ray_cast_dataset_save.bin";
  "QmXDu8XKrokoeGZg3wkL7r6zLMNyaiaHhS1cC3Fjjwkqoi", "incube_ray_cast_1m.bin";
  "QmQmp6usmpequDCg1Q1RBLNRUipJjoY5be9AZpa8e9YRdj", "incube_ray_cast_2m.bin";
  "QmcWMCqSZRu7Sg8mU8qsSMQNVq2ZEQLZukg1rpRRNVLxwS", "incube_ray_cast_5m.bin";
  "Qmeh5UKMjQwRvCqdYo92GVpSXe7FY1XkCieUQZj9K3C4bw", "incube_ray_cast_5m.txt";
  "QmR2FHtjt2kyYp1fUg6UbcjkFCXCNfqRwjt5evCUjPuMzo", "livejournal.bin";
  "QmQ5Hg634fcLcGnbzbEDBbSb89TByhpRDr65LK1V4T45Ee", "loop_1010_10.bin";
  "QmWFRhiofLxKpeoKwerFFuKtYLQnNE2g5ybEcbPnrjh89x", "loop_1010_1000.bin";
  "Qma8Wb4hq3EnzbeAJoNMs9JHpznxRbNBDjScLfdRZqDKsX", "loop_1010_100000.bin";
  "QmTswjsAsBmqPCubzNU2hBw5tC2JRDg8q11SDC9b3XCF4k", "loop_1010_1000000000.bin";
  "QmUzeR8yJbC6RrdKoyvnsaL49xVMXm51ZZrsnmMUAr8yqv", "loop_1010_33222591.bin";
  "QmethjsWxJC2aMvhavjgjmoivVm79v18jWkF1B11sceUJY", "loop_109_10.bin";
  "QmW1CT3D88B4uHkafKYPengNK93fa9xPUv7oMiVdqfGsd1", "loop_109_1000.bin";
  "QmcdfVa8YZQeqvYHxqrbCeypVYmQ73cEgyPpm9iTaPzWTD", "loop_109_100000000.bin";
  "QmYuxJVfgdG4cQqWBmXneanwj2WPMQTTDzZKcheAESNhwp", "loop_109_30000.bin";
  "Qmf4ds7ZRFhz3JrpvyM78rzqiAYJh9ecZ1v2rYs8v3GNVC", "loop_109_3000000.bin";
  "QmZkqdPc7qVVY3wZgYfA8XgSSjfHkL5xNhpA3Xd7gpkFBX", "mat27_large.bin";
  "QmZvyM8j6xRBZbwtKdzWYpR634QYkBnsQ8HsDK5Y5EHY4n", "mat27_medium.bin";
  "QmTtq7Fo2k8NNHSgYxLRtgzcNMnz5wgg8krFxWRLyHUTFG", "onsphere_ray_cast_1m.bin";
  "QmZYkpXs6eCwje6aDxn7d3K1DD13WaLggqWqCQDuNtQXkR", "onsphere_ray_cast_2m.bin";
  "QmQaJeP94yqZsQjLPgTJhNG4gNqprg3ZFNmYqpzFnqHfau", "onsphere_ray_cast_5m.bin";
  "QmdD3uK6hUfnK7W2pMpfFjgLfbJJYK6a95pXSe8ouyXVFb", "paths_100_phases_1_large.bin";
  "QmQYruSfm28CETx1WC2KX8gK2arFUrgkwaNrxHGFGQnwYD", "phased_524288_single_large.bin";
  "QmZnCAWbXDK9bsJN6Vgsad55kms8RGgH652Zg2mv64gV3J", "phased_low_50_large.bin";
  "QmUkCiAYaj3tufQ3B3WDPNYFrc77m7uMPaTUue9z3RrGL7", "phased_mix_10_large.bin";
  "QmfETQkR9PU2mwdFqis4qHHqfkb5AGjZjATK2TP5qyzncU", "rMatGraph_J_5_10000000.bin";
  "QmZX9N8iG4NXDAsBCSHhKtsBAonrisHsmPxnGX6Jkwp2B2", "randLocalGraph_J_5_10000000.bin";
  "QmNNunrwrWf1Ebz4xDAriQkNKpgmQbGXX1d2GG9vLX2neE", "random_arity_100_large.bin";
  "Qme1WH2x4qWzk7P2ojxpawnyAS15993M8jdh9YMrDSnxha", "refine_triangles_point2d_delaunay_in_square_large.bin";
  "QmckMba6Jkv7WjoULkweqhR4KQxMGatpsuYM5rg6vVU8js", "refinetriangles_point2d_delaunay_kuzmin_large.bin";
  "QmYGqukAsFSQxe8FEff9h65nGhCH23uSWwJKpH7YkCAijk", "rgg.bin";
  "QmQvdBNowHoHCn5LXRtxLs3aaEs8GLuz9yZk3HzF6jBA6x", "rmat24_large.bin";
  "QmUmE6UvxnxwNYAB1sgRx8RtRAzMYpBkbfagxomoJtNHZy", "rmat24_medium.bin";
  "QmesvvU9bKJkHpKx93F47PLvuVGnjs7M9mqj6sW3ZKsJB2", "rmat24_small.bin";
  "QmZv6vkPYwdoHimpDXBB9WrhomiyQ8asq1FyPKnwyzX53A", "rmat27_large.bin";
  "QmPNPzj9jTifjUntnLi71LoPCfmPWMWgZwJg6rgW8fLutX", "rmat27_medium.bin";
  "Qmdn85duXK1GsiQMRjQcdCaDopNqvLkitnY3T9abrkFNzE", "rmat27_small.bin";
  "QmUgS3XGE1C6AQJry1Vemxidocx7Bbc91sAsXygFneb4qN", "string_trigrams_large.bin";
  "QmPS8uF32prKEfc8hSUupnvj37fEiVTPqjzzwSZkeu3ny7", "string_trigrams_medium.bin";
  "QmVHyrpKstu7TWuYipftYrCCjExUJWncc6Gn9vHGKKN3iv", "string_trigrams_small.bin";
  "QmbNCimmwsy65Bg9ppFZyx5aDgCRiRNPJXMHhtHPypdetq", "tree_2_512_1024_large.bin";
  "QmX7h4VrE72Vb73dgHreaP2LFEbQxcXtxqovoL7WaJi6TG", "triangles_point2d_delaunay_in_square_medium.bin";
  "Qme1WH2x4qWzk7P2ojxpawnyAS15993M8jdh9YMrDSnxha", "triangles_point2d_delaunay_in_square_refine_large.bin";
  "Qme1WH2x4qWzk7P2ojxpawnyAS15993M8jdh9YMrDSnxha", "triangles_point2d_delaunay_in_square_small.bin";
  "QmUGcvjmXYyCxjhydHQzDc8uu6X5X4Uqjdvs4QFjfeMFMi", "triangles_point2d_delaunay_kuzmin_medium.bin";
  "QmckMba6Jkv7WjoULkweqhR4KQxMGatpsuYM5rg6vVU8js", "triangles_point2d_delaunay_kuzmin_refine_large.bin";
  "QmckMba6Jkv7WjoULkweqhR4KQxMGatpsuYM5rg6vVU8js", "triangles_point2d_delaunay_kuzmin_small.bin";
  "QmaRTv7FZdAi63nDKDxUV6goQLxRdpsM2H9dDkERPoAxfH", "turbine_ray_cast_dataset.bin";
  "QmcuLdFKKydtFvJWugg4v8AVBwsMcCMrSm5zKUjeWLcEaf", "twitter.bin";
  "QmRaZvZkbYkj67DEh5iXUQft3TX2RFy4G4GgD8VCqFG4cF", "unbalanced_tree_trunk_first_large.bin";
  "QmPfv38PWzfD9NKPsmTGRns8cQsYUgu3yTor4xeSbi4Vgt", "usa.bin";
  "Qmepjnfng9sMYbGrZt63MRZcS7g2WdB5D9pwhRdDK83etJ", "wikipedia-20070206.bin";
  "QmfYr68dj2CPKvf7Rg44Ae1yKXpUVTgyv2MJb3yhsX3UE2", "wikisamp.xml.bin";
  "QmWNiBmwPn7herCEjnc8b7iz6GZRneYGuoQWmsZEH4jU4Z", "xyzrgb_dragon_ray_cast_dataset.bin";
  "QmT76cbeEP64rS617yXpXr3efaAtBXa8mUZ4PTBVNuUYMd", "xyzrgb_manuscript_ray_cast_dataset.bin";
]

let row_of_infile path_to_infile infile =
  let h, _ = List.find (fun (_, f) -> (path_to_infile f) = infile) infiles_by_hash in
  (h, infile)

let infile_of_input_descriptor (p, _, _) = p

let fetch_infiles_of path_to_infile force_get is_virtual descrs =
  let infiles = List.map infile_of_input_descriptor descrs in
  let table = List.map (row_of_infile path_to_infile) infiles in
  ipfs_get_files table force_get is_virtual

let fetch_infiles_of_benchmark path_to_infile force_get is_virtual (benchmark : benchmark_descriptor) =
  fetch_infiles_of path_to_infile force_get is_virtual benchmark.bd_input_descr

let fetch_infiles_of_benchmarks path_to_infile force_get is_virtual all_benchmarks benchmarks =
  let keep_benchmark (benchmark : benchmark_descriptor) = List.exists (fun n -> n = benchmark.bd_name) benchmarks in
  let selected_benchmarks = List.filter keep_benchmark all_benchmarks in
  List.iter (fetch_infiles_of_benchmark path_to_infile force_get is_virtual) selected_benchmarks

let benchmarks =
  let p b =
    List.exists (fun a -> b.bd_name = a) all_benchmarks
  in
  List.filter p benchmarks'

let input_descriptors =
  List.flatten (List.map (fun b -> b.bd_input_descr) benchmarks)

let pretty_input_name n =
  match List.find_all (fun (m, _, _) -> m = n) input_descriptors with
  | (m, _, p) :: _ -> p
  | [] -> failwith ("pretty name: " ^ n)

let get() = 
  let _ = system (Printf.sprintf "mkdir -p %s" arg_path_to_data) arg_virtual_get in    
  (
    fetch_infiles_of_benchmarks path_to_infile arg_force_get arg_virtual_get benchmarks all_benchmarks;
    ())
                  
let make() =
  build "." all_progs arg_virtual_build

let mk_never_promote =
  mk int "never_promote" 1

let file_results_heartbeat_elision exp_name =
  file_results (exp_name ^ "_heartbeat_elision")

let file_results_pbbs_elision exp_name =
  file_results (exp_name ^ "_cilk_elision")

let file_results_heartbeat_single_proc exp_name =
  file_results (exp_name ^ "_heartbeat_single_proc")

let file_results_pbbs_single_proc exp_name =
  file_results (exp_name ^ "_pbbs_single_proc")

let nb_proc = List.length arg_proc
let nb_multi_proc = List.length multi_proc
        
let run() =
  List.iter (fun benchmark ->
    let r mk_progs file_results = 
      Mk_runs.(call (run_modes @ [
        Output file_results;
        Timeout 400;
        Args (mk_progs & benchmark.bd_infiles); ]))
    in
    let heartbeat_prog = mk_heartbeat_prog benchmark.bd_name in
    let heartbeat_elision_prog = heartbeat_prog & mk_never_promote in
    let pbbs_prog = mk_pbbs_prog benchmark.bd_name in
    let pbbs_elision_prog = mk_pbbs_elision_prog benchmark.bd_name in
    (if nb_multi_proc > 0 then (
      r ((heartbeat_prog ++ pbbs_prog) & mk_multi_proc) (file_results benchmark.bd_name))
     else
       ());
    (if List.exists (fun p -> p = 1) arg_proc then (
      r (heartbeat_prog & mk_single_proc) (file_results_heartbeat_single_proc benchmark.bd_name);
      r (pbbs_prog & mk_single_proc) (file_results_pbbs_single_proc benchmark.bd_name);
      r (heartbeat_elision_prog & mk_single_proc) (file_results_heartbeat_elision benchmark.bd_name);
      r (pbbs_elision_prog & mk_single_proc) (file_results_pbbs_elision benchmark.bd_name))
     else
       ())
  ) benchmarks

let check = nothing  (* do something here *)

let plot() =
    let tex_file = file_tables_src name in
    let pdf_file = file_tables name in

    let main_formatter =
      Env.format (Env.(
                  [
                   ("proc", Format_custom (fun n -> ""));
                   ("lib_type", Format_custom (fun n -> ""));
                   ("infile", Format_custom pretty_input_name);
                   ("prog", Format_custom (fun n -> ""));
                   ("type", Format_custom (fun n -> ""));
                   ("source", Format_custom (fun n -> ""));
                 ]
                 ))
    in
    let nb_application_cols = 2 in
    let nb_seq_elision_cols = 2 in
    let nb_single_core_cols = 2 in
    let nb_multi_core_cols = 4 + (if arg_show_utilization then 2 else 0) in
    let nb_cols = nb_application_cols + nb_seq_elision_cols + nb_single_core_cols + (nb_multi_proc * nb_multi_core_cols) in

    Mk_table.build_table tex_file pdf_file (fun add ->
      let hdr =
        let ls = String.concat "|" (XList.init (nb_cols - 1) (fun _ -> "c")) in
        Printf.sprintf "|p{1cm}l|%s" ls
      in
      add (Latex.tabular_begin hdr);

      (* Emit first row, i.e., first-level column labels *)
      Mk_table.cell ~escape:true ~last:false add (Latex.tabular_multicol nb_application_cols "|l|" "Application/input");
      Mk_table.cell ~escape:true ~last:false add (Latex.tabular_multicol nb_seq_elision_cols "|l|" "Sequential elision");
      Mk_table.cell ~escape:true ~last:false add (Latex.tabular_multicol nb_single_core_cols "|c|" "1-core execution");
      ~~ List.iteri multi_proc (fun i proc ->
        let last = i + 1 = nb_multi_proc in
	      let label = Printf.sprintf "%d-core execution" proc in
        Mk_table.cell ~escape:false ~last:last add (Latex.tabular_multicol nb_multi_core_cols "c|" label));
      add Latex.tabular_newline;

      (* Emit second row, i.e., second-level column labels *)
      for i = 1 to nb_application_cols do
        Mk_table.cell ~escape:false ~last:false add ""
      done;
      Mk_table.cell ~escape:false ~last:false add "PBBS";
      Mk_table.cell ~escape:false ~last:false add "Encore";
      Mk_table.cell ~escape:false ~last:false add "PBBS";
      Mk_table.cell ~escape:false ~last:false add "Encore";
      ~~ List.iteri multi_proc (fun i proc ->
        let last = i + 1 = nb_multi_proc in
	      Mk_table.cell ~escape:false ~last:false add "PBBS";
	      Mk_table.cell ~escape:false ~last:false add "Encore";
	      if arg_show_utilization then begin
  	        Mk_table.cell ~escape:false ~last:false add "PBBS";
		Mk_table.cell ~escape:false ~last:false add "Encore"
              end;
	      Mk_table.cell ~escape:false ~last:last add (Latex.tabular_multicol 2 "|c|" "Encore/PBBS"));
      add Latex.tabular_newline;

      (* Emit third row, i.e., third-level column labels *)
      for i = 1 to nb_application_cols do
        Mk_table.cell ~escape:false ~last:false add ""
      done;
      Mk_table.cell ~escape:false ~last:false add "(s)";
      Mk_table.cell ~escape:false ~last:false add "";
      Mk_table.cell add (Latex.tabular_multicol 2 "|l|" "(relative to elision)");
      ~~ List.iteri multi_proc (fun i proc ->
        let last = i + 1 = nb_multi_proc in
	      Mk_table.cell ~escape:false ~last:false add "(s)";
	      Mk_table.cell ~escape:false ~last:false add "";
	      if arg_show_utilization then begin
	        Mk_table.cell ~escape:false ~last:false add (Latex.tabular_multicol 2 "|c|" "Utilization")
              end;
	      Mk_table.cell ~escape:false ~last:false add "Idle time";
	      Mk_table.cell ~escape:false ~last:last add "Nb. threads");
      add Latex.tabular_newline;

      (* Emit two rows for each benchmark *)
      ~~ List.iteri benchmarks (fun benchmark_i benchmark ->
        Mk_table.cell add (Latex.tabular_multicol nb_application_cols "|l|" (sprintf "\\textbf{%s}" (Latex.escape benchmark.bd_name)));
	      let nbc = nb_cols - nb_application_cols in
        for i = 1 to nbc do
          let last = i = nbc in
          Mk_table.cell ~escape:true ~last:last add "";
        done;
        add Latex.tabular_newline;
        let results_file_pbbs_elision = file_results_pbbs_elision benchmark.bd_name in
        let results_pbbs_elision = Results.from_file results_file_pbbs_elision in
        let results_file_heartbeat_elision = file_results_heartbeat_elision benchmark.bd_name in
        let results_heartbeat_elision = Results.from_file results_file_heartbeat_elision in
        let results_file_pbbs_single_proc = file_results_pbbs_single_proc benchmark.bd_name in
        let results_pbbs_single_proc = Results.from_file results_file_pbbs_single_proc in
        let results_file_heartbeat_single_proc = file_results_heartbeat_single_proc benchmark.bd_name in
        let results_heartbeat_single_proc = Results.from_file results_file_heartbeat_single_proc in
	      let results_file = file_results benchmark.bd_name in
	      let all_results = Results.from_file results_file in
	      let results = all_results in
	      let env = Env.empty in
	      let env_rows = benchmark.bd_infiles env in
        ~~ List.iter env_rows (fun env_rows ->  (* loop over each input for current benchmark *)
          let results = Results.filter env_rows results in
          let results_pbbs_single_proc = Results.filter env_rows results_pbbs_single_proc in
          let results_heartbeat_single_proc = Results.filter env_rows results_heartbeat_single_proc in
          let env = Env.append env env_rows in
          let input_name = main_formatter env_rows in
          let _ = Mk_table.cell ~escape:true ~last:false add "" in
          let _ = Mk_table.cell ~escape:true ~last:false add input_name in
	        let pbbs_elision_sec =
            let results_pbbs_elision = Results.filter env_rows results_pbbs_elision in
	          let [col] = ((mk_pbbs_elision_prog benchmark.bd_name) & mk_single_proc) env in
	          let results = Results.filter col results_pbbs_elision in
	          Results.get_mean_of "exectime" results
	        in
	        let heartbeat_elision_sec =
            let results_heartbeat_elision = Results.filter env_rows results_heartbeat_elision in
            let [col] = (mk_heartbeat_prog benchmark.bd_name & mk_never_promote & mk_single_proc) env in
	          let results = Results.filter col results_heartbeat_elision in
	          Results.get_mean_of "exectime" results
	        in
	        let heartbeat_elision_rel_pbbs_elision = string_of_percentage_change pbbs_elision_sec heartbeat_elision_sec in
	        let _ = (
            Mk_table.cell ~escape:false ~last:false add (Printf.sprintf "%.2f" pbbs_elision_sec);
            Mk_table.cell ~escape:false ~last:false add heartbeat_elision_rel_pbbs_elision)
	        in
	        let pbbs_single_proc_sec =
            let results_pbbs = Results.filter env_rows results_pbbs_single_proc in
	          let [col] = ((mk_pbbs_prog benchmark.bd_name) & mk_single_proc) env in
	          let results = Results.filter col results_pbbs in
	          Results.get_mean_of "exectime" results
	        in
	  let heartbeat_single_proc_sec =
      let results_heartbeat = Results.filter env_rows results_heartbeat_single_proc in
	    let [col] = ((mk_heartbeat_prog benchmark.bd_name) & mk_single_proc) env in
	    let results = Results.filter col results_heartbeat in
	    Results.get_mean_of "exectime" results
	  in
	  let pbbs_single_proc_rel_pbbs_elision = string_of_percentage_change pbbs_elision_sec pbbs_single_proc_sec in
	  let heartbeat_single_proc_rel_heartbeat_elision = string_of_percentage_change heartbeat_elision_sec heartbeat_single_proc_sec in
	  let _ = (
      Mk_table.cell ~escape:false ~last:false add pbbs_single_proc_rel_pbbs_elision;
      Mk_table.cell ~escape:false ~last:false add heartbeat_single_proc_rel_heartbeat_elision)
	  in
    ~~ List.iteri multi_proc (fun proc_i proc ->
      let last = proc_i + 1 = nb_multi_proc in
      let mk_procs = mk int "proc" proc in
	    let (pbbs_sec, pbbs_utilization, pbbs_idle_time, pbbs_multi_proc_nb_threads) =
        let [col] = ((mk_pbbs_prog benchmark.bd_name) & mk_procs) env in
        let env = Env.append env col in
        let results = Results.filter col results in
        let sec = eval_exectime env all_results results in
	      let util = Results.get_mean_of "utilization" results in
	      let idle_time = util *. sec in
	      let nb_threads = Results.get_mean_of "nb_threads_alloc" results in
	      let nb_threads = if nb_threads = 0. then 1. else nb_threads
	      in
  	    (sec, util, idle_time, nb_threads)
      in
	    let (heartbeat_sec, heartbeat_utilization, heartbeat_idle_time, heartbeat_multi_proc_nb_threads) =
        let [col] = ((mk_heartbeat_prog benchmark.bd_name) & mk_procs) env in
        let env = Env.append env col in
        let results = Results.filter col results in
        let sec = eval_exectime env all_results results in
	      let util = Results.get_mean_of "utilization" results in
	      let idle_time = util *. sec in
 	      let nb_threads = Results.get_mean_of "nb_promotions" results in
	      let nb_threads = if nb_threads = 0. then 1. else nb_threads
	      in
  	    (sec, util, idle_time, nb_threads)
      in
	    let heartbeat_rel_pbbs = string_of_percentage_change pbbs_sec heartbeat_sec in
	    let pbbs_utilization_str = string_of_percentage ~show_plus:false pbbs_utilization in
	    let heartbeat_utilization_str = string_of_percentage ~show_plus:false heartbeat_utilization in
	    let idle_time_enc_by_pbbs_str = string_of_percentage_change pbbs_idle_time heartbeat_idle_time in
	    let nb_threads_enc_by_pbbs_str = string_of_percentage_change pbbs_multi_proc_nb_threads heartbeat_multi_proc_nb_threads in
	    Mk_table.cell ~escape:false ~last:false add (Printf.sprintf "%.2f" pbbs_sec);
	    Mk_table.cell ~escape:false ~last:false add heartbeat_rel_pbbs;
	    if arg_show_utilization then begin
	      Mk_table.cell ~escape:false ~last:false add pbbs_utilization_str;
	      Mk_table.cell ~escape:false ~last:false add heartbeat_utilization_str
	    end;
	    Mk_table.cell ~escape:false ~last:false add idle_time_enc_by_pbbs_str;
	    Mk_table.cell ~escape:false ~last:last add nb_threads_enc_by_pbbs_str);
    add Latex.tabular_newline);
  );
  add Latex.tabular_end;
  add Latex.new_page;
  ())

let all () = select get make run check plot

end
    
(*****************************************************************************)
(** Main *)

let _ =
  let arg_actions = XCmd.get_others() in
  let bindings = [
    "sequence", ExpSequenceLibrary.all;
    "compare", ExpCompare.all;
  ]
  in
  Pbench.execute_from_only_skip arg_actions [] bindings;
  ()
