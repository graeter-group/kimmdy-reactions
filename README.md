# A collection of simple reactions for KIMMDY

[![test latest release](https://github.com/graeter-group/kimmdy-reactions/actions/workflows/tests.yml/badge.svg?branch=release-please--branches--main)](https://github.com/graeter-group/kimmdy-reactions/actions/workflows/tests.yml/?branch=release-please--branches--main)

## Documentation

See [KIMMDY documentation](https://graeter-group.github.io/kimmdy/)
and [documentation for these reactions](https://graeter-group.github.io/kimmdy-reactions/).

## Installation

Together with KIMMDY

```bash
pip installl kimmdy[reactions]
```

To install it separatly:

```bash
pip install kimmdy-reactions
```

## Making your own

- Implement your reaction as a subclass of `kimmdy.reaction.Reaction`
- Register your Reaction class in the **[options.entry_points]** section in the setup.cfg.
  The name you give here must match the entry in the config.yml for Kimmdy!

