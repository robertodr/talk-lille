with import (builtins.fetchGit {
  name = "nixos-19.03-2019-06-26";
  url = "https://github.com/NixOS/nixpkgs-channels";
  ref = "nixos-19.03";
  # Commit hash for nixos-19.03 as of 2019-06-26
  # `git ls-remote https://github.com/nixos/nixpkgs-channels nixos-19.03`
  rev = "75a88c1b9d0abbee53ce0271a4c360219a99787e";
}) {
  overlays = [(self: super: {})];
};

stdenv.mkDerivation {
  name = "talk-lille";
  buildInputs = [
    pipenv
  ];

  src = null;
  shellHook = ''
  SOURCE_DATE_EPOCH=$(date +%s)
  '';
}
