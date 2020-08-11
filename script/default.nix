{ pkgs   ? import <nixpkgs> {},
  stdenv ? pkgs.stdenv,
  sources ? import ./default-sources.nix,
  pbench ? import ./bench-script.nix { sources = sources; },
  cmdline ? sources.cmdline,
  chunkedseq ? sources.chunkedseq,
  cilk-plus-rts-with-stats ? sources.cilk-plus-rts-with-stats,
  pbbs-include ? sources.pbbs-include,
  cactus-stack ? sources.cactus-stack,
  sptl ? sources.sptl,
  heartbeatSrc ? sources.heartbeatSrc,
  gperftools ? pkgs.gperftools,
  hwloc ? pkgs.hwloc,
  libunwind ? null,
  gcc ? pkgs.gcc7,
  clang ? null,
  pathToResults ? "",
  pathToData ? "",
  nbBuildCores ? 0,
  buildDocs ? false
}:

stdenv.mkDerivation rec {
  name = "heartbeat-${version}";
  version = "v1.0-pldi";

  src = heartbeatSrc;

  buildInputs =
    let docs = if buildDocs then [ pkgs.pandoc ] else []; in
    [ sptl pbbs-include cmdline chunkedseq
      pkgs.makeWrapper pkgs.R pkgs.texlive.combined.scheme-small
      pkgs.ocaml gcc hwloc
    ] ++ docs ++ (if clang != null then [ clang ] else []);
        
  buildPhase =
    let hwlocConfig =
      if hwloc != null then 
        ''USE_HWLOC=1 \
          USE_MANUAL_HWLOC_PATH=1 \
          HWLOC_FLAGS="-I ${hwloc.dev}/include/" \
          HWLOC_LIBS="-L ${hwloc.lib}/lib/ -lhwloc"
      '' else "";
    in
    let clangConfig =
          if clang != null then ''
          export CXX=${clang}/clang++ -frtti -fno-exceptions -latomic
          '' else "";
    in
    let sptlConfigFile = pkgs.writeText "sptl_config.txt" "${sptl}/bin/"; in
    ''
    cp ${sptlConfigFile} bench/sptl_config.txt
    ${clangConfig}
    make -C bench clean
    make -j $NIX_BUILD_CORES -C bench \
      all \
      CHUNKEDSEQ_PATH=${chunkedseq}/include/ \
      CACTUS_PATH=${cactus-stack}/include/ \
      CMDLINE_PATH=${cmdline}/include/ \
      SPTL_PATH=${sptl}/include/ \
      PBBS_INCLUDE_PATH=${pbbs-include}/include/ \
      HEARTBEAT_INCLUDE_PATH=../include/ \
      USE_CILK=1 \
      CUSTOM_MALLOC_PREFIX="-ltcmalloc -L${gperftools}/lib" \
      CILK_EXTRAS_PREFIX="-L ${cilk-plus-rts-with-stats}/lib -I ${cilk-plus-rts-with-stats}/include -ldl -DCILK_RUNTIME_WITH_STATS -DCUSTOM_CILK_PLUS_RUNTIME" \
      ${hwlocConfig}
  '';
  
  installPhase =
    let lu =
        if libunwind != null then
           ''--prefix LD_LIBRARY_PATH ":" ${libunwind}/lib''
        else "";
    in
    let hw =
        if hwloc != null then
          ''--prefix LD_LIBRARY_PATH ":" ${hwloc.lib}/lib''
        else "";
    in
    let nmf = "-skip make";
    in
    let rf =
      if pathToResults != "" then
        "-path_to_results ${pathToResults}"
      else "";
    in
    let df =
      if pathToData != "" then
        "-path_to_data ${pathToData}"
      else "";
    in
    let flags = "${nmf} ${rf} ${df}";
    in
    ''
    mkdir -p $out/bench/
    cp ${pbench}/bench $out/bench/bench.pbench
    wrapProgram $out/bench/bench.pbench --prefix PATH ":" ${pkgs.R}/bin \
       --prefix PATH ":" ${pkgs.texlive.combined.scheme-small}/bin \
       --prefix PATH ":" ${gcc}/bin \
       --prefix PATH ":" ${pkgs.ipget}/bin \
       --prefix PATH ":" $out/bench \
       --prefix LD_LIBRARY_PATH ":" ${gcc}/lib \
       --prefix LD_LIBRARY_PATH ":" ${gcc}/lib64 \
       --prefix LD_LIBRARY_PATH ":" ${gperftools}/lib \
       --prefix LD_LIBRARY_PATH ":" ${cilk-plus-rts-with-stats}/lib \
       --set TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD 100000000000 \
       ${lu} \
       ${hw} \
       --add-flags "${flags}"
    make -C bench INSTALL_FOLDER=$out/bench install
    cp bench/sptl_config.txt $out/bench/sptl_config.txt
    cp bench/nb_cores $out/bench/
    '';

  meta = {
    description = "The C++ library for benchmarking the Heartbeat Scheduling algorithm for fine-grain, irregular parallelism.";
    license = "MIT";
    homepage = http://deepsea.inria.fr/heartbeat/;
  };
}
