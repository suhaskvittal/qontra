# QontraSim ISA

This file summarizes all valid QontraSim instructions that can be used in QES
and the `FullSystemSimulator` architecture.

## Gates

Valid gates: `h`, `cx`, `s`, `x`, `z`, `sdg`, `measure`, `reset`. Gates are
specified as SIMD-like instructions, much like how Stim specifies gates. For
instance, `h 0,1,2,3` specifies a **H** gate on qubits 0 through 3. Similarly,
`cx 0,1,2,3` specifies **CX**(0, 1) and **CX**(2, 3).

## Events and observables

Detection events can be specified via the syntax:
```
event <event-no>, <meas-no>, <meas-no>, ...;
```
The corresponding detection event is written to syndrome bit `<event-no>` and
the bit is the XOR of the measurement bits specified. Observables are similar,
except the first argument specifies the observable to write to. In
`full_system_sim.h`, `syndrome_table[*][event-no]` and
`observable_table[*][obs-no]` are set.

## Tracking Management Instructions

Much like a real system, the `FullSystemSimulator` has a limited amount of
record space for measurements, syndromes, and observables. Consequently, reusing
this space once information is no longer needed can be useful. By default, 8
Kbits of storage per trial is allocated for all tracking structures. This value
can be updated by setting `G_RECORD_SPACE_SIZE` prior to making a simulator
object.

Two instructions are supported to make additional space. As older data is
generally stale and unneeded, both instructions are shifts:
```
mshift <by>;
eshift <by>;
```
For `mshift` (`eshift`), the measurement (event) counter is shifted to the left
(reduced) by the operand. For instance, if the measurement counter is at 25 and
`mshift 4` is executed, the new counter will be at 21.

Shifting is one method of tracking management. We also provide offset instructions:
```
moffset <offset>;
eoffset <offset>;
```
While `mshift` and `eshift` update the recorded values, `moffset` and `eoffset`
change where data is being recorded. For instance, if an event would normally be
recorded at position *X*, after calling `eshift 4`, it would now be recorded at
position *X*+4. Note that `moffset` and `eoffset` are cumulative, so an additional
call `eshift 3` would result in the event being recorded at *X*+7, not *X*+3.
Note that the offset change can be undone by providing a negative offset.

## Branching Instructions

`FullSystemSimulator` implements two types of branching instructions:
1. Permissive branches, which can be used with a batch size greater than one.
2. Prohibitive branches, which cannot be used with a batch size greater than
   one.
The reason for this decision is to avoid code bloat, where one simulator is
exclusively used for batch simulation, and another is used for trial-by-trial
simulation. Having two simulators hinders maintainability, especially as the
functions and data structures in both simulators would likely be the same.

Nevertheless, we detail what each type of branch entails.

### Permissive Branches

Permissive branches are generally compatible with batching and do not require
gymnastics to implement. Such instructions are:
```
branch_and_wait <address>, <register>;
return_if <register>;
```

`branch_and_wait` is an instruction that is a blocking branch. If the specified
register contains a nonzero value (see [here](#registers) for a list of
registers), then the PC for the trial jumps to the address. Then, the simulator
waits until all trials reach the same address before proceeding. Obviously, such
an implementation of the instruction is *not* perfect, as a program can have a
deadlock with the following instructions:
```
...
branch_and_wait     DEAD, $r1;   // Trial A branches here
...
branch_and_wait     BEEF, $r2;   // Trial B branches here
...
BEEF: nop
DEAD: nop
```
As trial A will never reach `BEEF` since it is stuck waiting at `DEAD`, the
simulator stalls at `BEEF`. Thus, we implement `branch_and_wait` to continue if
all trials have their PC *greater than or equal* to the target address. Thus, in
the above example, the deadlock is avoided as all trials (except trial A) will
reach `BEEF`, and then they will jointly proceed to `DEAD`, at which point trial
A can proceed with the computation. This is suitable for most use cases in QEC,
as the PC generally only increases.

`return_if` is a special instruction for subroutines, where if the register's
contents are nonzero, the microcode is exited. If `return_if` succeeds, then
the condition register `$rbrk` is set.

### Prohibitive Branches

## Bitwise Operations

The only operations supported on registers are bitwise operations. These are:
```
not <register>;                 // A = ~A
or  <register>, <register>;     // A = A | B
and <register>, <register>;     // A = A & B
xor <register>, <register>;     // A = A ^ B
```

## Data Movement

Memory in the `FullSystemSimulator` is split between that of the register file
and the dedicated memory for measurements, syndromes, and observables. To
facilitate data movement between registers and between registers and dedicated
memory, we provide the following instructions:
```
mov     <register>, <register>;
movm    <register>, <meas-no>; 
move    <register>, <event-no>;
movo    <register>, <obs-no>;
```

## Registers

General purpose registers are `$r1` through `$r32`, and `$r0` is a special zero
register. Special registers are as follows:
1. `$rbrk` which is a condition register for the `return_if` instruction.

As general-purpose arithmetic is not supported in `FullSystemSimulator`, all
registers are bitvectors currently and simply store one bit.

## Other

A program can call the execution of a subroutine with the instruction
```
call <subroutine>;
```
The program then steps into the code of the subroutine and returns to the
calling function once the subroutine has finished. Such a subroutine may be a
syndrome extraction function, for example.
