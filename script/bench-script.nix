{ pkgs   ? import <nixpkgs> {},
  stdenv ? pkgs.stdenv,
  makeWrapper ? pkgs.makeWrapper,
  hwloc ? pkgs.hwloc,
  ipget ? pkgs.ipget,
  sources ? import ./default-sources.nix,
  buildDunePackage ? pkgs.ocamlPackages.buildDunePackage,
  pbenchOcaml ? import sources.pbenchOcamlSrcs.pbenchOcaml { pbenchOcamlSrc = sources.pbenchOcamlSrcs.pbenchOcamlSrc; }
}:

let benchDune =
      buildDunePackage rec {
        pname = "bench";
        version = "1.0";
        src = sources.benchOcamlSrc;
        buildInputs = [ pbenchOcaml ];
      };
in

let bench = import sources.pbenchOcamlSrcs.pbenchCustom {
                      benchSrc = sources.benchOcamlSrc;
                      bench = "${benchDune}/bin/bench"; };
in

stdenv.mkDerivation rec {
  name = "bench-script";

  src = "${sources.heartbeatSrc}";

  buildInputs = [ makeWrapper ];

  configurePhase =
    let getNbCoresScript = pkgs.writeScript "get-nb-cores.sh" ''
      #!/usr/bin/env bash
      nb_cores=$( ${hwloc}/bin/hwloc-ls --only core | wc -l )
      echo $nb_cores
    '';
    in
    ''
    cp ${getNbCoresScript} get-nb-cores.sh
    '';

  installPhase = ''
    mkdir -p $out
    cp get-nb-cores.sh $out/
    cp ${bench}/bench $out/bench
    wrapProgram $out/bench \
      --prefix PATH ":" $out/ \
      --prefix PATH ":" ${ipget}/bin
    '';


}
