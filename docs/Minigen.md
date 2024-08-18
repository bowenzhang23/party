# Minigen

## Overview

Write a minimal parton-level event generator in Rust.

It can do ee -> y -> mumu at a lepton collider and qqbar -> Z/y -> mumu and a hadron collider.

It also contains **physics constants** used in particle physics

Based on the reference: https://arxiv.org/pdf/1412.4677.pdf

> This is not a general purpose event generator

## Devlog

### Basic

[done] ee -> y -> mumu
[done] qqbar -> Z/y -> mumu (q-q interaction)

### Advanced

[done] parton distribution function
[done] Support multi-D integration
[done] qqbar -> Z/y -> mumu, (p-p interaction)
