Common L1T utilities for demonstrators
======================================

This package contains a prototype implementation of subsystem-agnostic utility functions and
classes for converting between emulator input/output objects (i.e. C++ objects representing
tracks, calo clusters, muon TPs, vertices etc) and the I/O buffer text files (i.e. the files
that are used to store data that is loaded into / captured from FPGAs, and used to play/capture
data in HDL simulations). The motivation for these tools and their scope was briefly summarised
in the L1T meeting on 16th February (see slides [here](https://indico.cern.ch/event/1008519/contributions/4234188/attachments/2191061/3703176/l1tPhase2_cmsswBufferIO_20210216.pdf))

A brief tour of the code:

 * One example EDAnalyzer - `GTTInputFileWriter` - that creates text files for loading into
   input buffers of the GTT board
   * Source: `plugins/GTTInputFileWriter.cc`
   * Test cmsRun config: `test/gtt/createFirmwareInputFiles_cfg.py`
 * One example EDProducer - `GTTOutputFileReader` - that reads text files that would be 
   produced by GTT vertex-finding FW
   * Source: `plugins/GTTOutputFileReader.cc`
   * Test cmsRun config: `test/gtt/verifyFirmwareOutput_cfg.py`
 * Main utility classes:
    - `BoardData`: Represents the data stored in I/O buffers
    - `ChannelSpec`: Simple struct containing parameters that define a link's packet
      structure (e.g. TMUX period/index, gap between packets)
    - `BoardDataReader`: This class ...
        1. reads a set of buffer files
        2. verifies data conforms to expected structure (largely defined by `map<channelIndex, ChannelSpec>`)
        3. splits out each event; and
        4. returns data for each event in turn via `getNextEvent()` method
    - `BoardDataWriter`: Essentially, does the opposite of the reader -
        1. accepts per-event data via `addEvent` method
        2. concatenates events according to specified structure; and
        3. automatically writes out pending data to file whenever the limit on the
           number of frames per file is reached.

(Note: For simplicity this code has been put in its own package during development,
but this might not be its final location.)

Given the above, the contents of `GTTInputFileWriter.cc` and `GTTOutputFileReader.cc`
should be mostly self-explanatory, but a couple of additional notes:

 * These `.cc` files each contain some hardcoded constants defining the link TM periods,
   TM indices, packet structure etc. that are used to create the `BoardDataReader`/`BoardDataWriter`
   instance (so that it can correctly separate/concatenate events).
   * I think at least some (if not all) of these types of constants should eventually be read
     from config files, to avoid e.g. needing to recompile code when the latency of an algo changes
 * The EDAnalyzer/EDProducer analyze/produce methods use functions from the 'codecs' directory
   to convert between EDM collections and vectors of Frames (the latter being one of the main
   building blocks of the `BoardData` instances that are returned by a `BoardDataReader` or fed
   to a `BoardDataWriter`).

Known areas for improvement
---------------------------

It would be good if people writing EDAnalyzers and EDProducers for systems where the link TMUX
period is different from the board TMUX period wouldn't have to keep track of `eventIndex % tmux`
for each event in order to put the data on the correct channel index. For example, I've had to do
that in the `GTTInputFileWriter` since the track links are TMUX 18, but the GTT is TMUX 6, and so 
depending on the event being processed the tracks from a specific eta-phi sector will come on one
of three links. One solution would be to:
 * Enhance the channel specs that are given to the `BoardDataWriter` and `BoardDataReader`
   constructors, such that the `BoardDataWriter::addEvent` & `BoardDataReader::getNextEvent` methods
   could return/receive a map of 'logical' event-index-independent channel indices (e.g. conceptually
   'track/N' where N=0,1,...,17); and
 * update the implementation of `BoardDataWriter::addEvent` & `BoardDataReader::getNextEvent` to
   handle the task of re-mapping between these logical event-index-independent channels and the
   set of channels that should be used for any given event.

Suggestions (and implementations) for improvements most welcome, but please get in touch with 
me about these before writing large amounts of code, to avoid duplicating efforts or divergent
solutions to the same problem.
