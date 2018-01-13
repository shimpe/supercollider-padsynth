(
//sample creation using the PadSynth algorith described at http://wiki.linuxmusicians.com/doku.php?id=zynaddsubfx_manual#padsynth_algorithm
//based on code from Donald Craig http://new-supercollider-mailing-lists-forums-use-these.2681727.n2.nabble.com/Epic-Pads-td7487382.html#a7492701
s.waitForBoot({
	var table, re, im, tab, fftsize, result, pars, freqs, amps, bandwidth, partials, note;
	var samplerate = s.sampleRate;
	var sgroup;
	var fxgroup1,fxgroup2;
	var reverbbus,chorusbus;
	var reverb, chorus;
	var prepareSingleBuffer, prepareAllBuffers;
	var spectrum, xvals;

	s.freeAll;
	s.freeAllBuffers;
	s.sync;

	sgroup = Group.new;
	fxgroup1 = Group.after(sgroup);
	fxgroup2 = Group.after(fxgroup1);
	reverbbus = Bus.audio(s, 1);
	chorusbus = Bus.audio(s, 1);

	~buffers = [];

	SynthDef(\padfilterenv, {
		| out=0, referencefreq=130, amp=0.5, freq=440, buffer=1, gate=1, attack=0.1, decay=0.2, sustain=0.6, release=2, filtercutoff=10000, filterresonance=1.0, filterattack=0.01, filterdecay=0.2, filtersustain=0.8, filterrelease=2, filtergain=1.0, glissando=0 |
		var env = EnvGen.ar(Env.adsr(attack, decay, sustain, release), gate, doneAction:Done.freeSelf);
		var env2 = EnvGen.ar(Env.adsr(filterattack, filterdecay, filtersustain, filterrelease), gate, doneAction:Done.none);
		var frequency = VarLag.ar(freq, glissando, warp:\exponential);
		var sig = env*PlayBuf.ar(1, buffer, rate:((frequency.cpsmidi)-(referencefreq.cpsmidi)).midiratio, loop:1);
		sig = RLPF.ar(sig, env2*filtercutoff, rq:filterresonance, mul:filtergain);
		Out.ar(out, amp*sig);
	}, rates:[nil, nil, nil, \ar, nil, nil, nil, nil, nil, nil, nil]).add;

	SynthDef(\reverb, {
		| out=0, inbus, mix=0.5, room=1.0 |
		var insig = FreeVerb.ar(In.ar(inbus, 1),mix,room);
		Out.ar(out, insig!2);
	}).add;

	SynthDef(\chorus, {
		| outbus=0, inbus, predelay=0.08, speed=0.05, depth=0.1, ph_diff=0.5 |
		var in, sig, modulators, numDelays = 12;
		in = In.ar(inbus, 1) * numDelays.reciprocal;
		modulators = Array.fill(numDelays, { |i|
			LFPar.kr(speed * rrand(0.94, 1.06), ph_diff * i, depth, predelay);
		});
		sig = DelayC.ar(in, 0.5, modulators);
		sig = sig.sum; //Mix(sig);
		Out.ar(outbus, sig!2); // output in stereo
	}).add;

	s.sync;

	// calculate single wavetable
	prepareSingleBuffer = {
		| partials /* flat list of [partial idx, partial amplitude, partial idx, partial amplitude, ...]. Partials often are integers (or close to) */,
		  min_length, /* min length of generated wave table in seconds */
		  spread /* band width used to generate new partials around the existing partials */,
		  reference_note /* note for which this spectrum is being generated */|

		var fftsize = (min_length*s.sampleRate).nextPowerOfTwo;
		var pars = (partials.size/2);
		var bandwidth = (1+spread);
		var note = reference_note;
		var	table = Signal.newClear(fftsize);
		var tab = Signal.fftCosTable(fftsize);
		var re = Signal.newClear(fftsize);
		var im = Signal.newClear(fftsize);
		var freqs = Array.newClear(pars);
		var amps = Array.newClear(pars);
		var buffer;
		var deinterlaced;
		fftsize.do({ |i|
			re[i] = 0.0;
			im[i] = 0.0;
			table[i] = 0.0;
		});

		// partials are specified in a flat list containing
		// each time partial number followed by corresponding partial volume.
		// first deinterlace this flat list into a list of frequencies and a list of amplitudes
		deinterlaced = partials.unlace;
		freqs = deinterlaced[0]*(note.midicps);
		amps = deinterlaced[1];

		// next we're going to generate extra (smeared) partials. This helps in adding life and warmth to the sounds.
		// if you specify spread == 0, no extra partials will be added
		pars.do({ |i|
			var freq, lo, hi,amp;
			freq = freqs[i];
			amp = amps[i];
			lo = ((freq/bandwidth)*(fftsize/samplerate)).round; // partial at frequency freq will be smeared over frequencies lo to hi
			hi = ((freq*bandwidth)*(fftsize/samplerate)).round;
			// generate extra partials between frequencies lo = freq/(1+spread) and hi = freq*(1+spread)
			(hi-lo+1).do({ |j|
				var mag, phase, val;
				var index = j.linlin(0, hi-lo, lo, hi);
				// only fill up lower half of spectrum:
				// right half later is derived from this left half to ensure a real-valued inverse fourier transform
				if(index < (fftsize/2), {
					if ((hi == lo), {
						mag = amp;
						table[index] = table[index] + mag; // add it to the result table
					}, /* else */ {
						val = j.linlin(0, hi-lo, -1, 1);
						mag = exp(val*val*10.0.neg) * amp; // generates a bell-shaped curve for val in [-1,1] with y-values between [-amp, amp]
						table[index] = table[index] + mag; // add it to the result table to create a "smeared" partial
					});
					phase = rrand(-pi, pi); // set random phase
					re[index] = re[index] + (cos(phase)*mag);
					im[index] = im[index] + (sin(phase)*mag);
				});
			});
		});

		// at this point, table contains the sum of all the specified + extra generated partials

		// calculate right half of spectrum to get a real-valued inverse FFT
		// right half must be the mirrored complex conjugate (i.e. make imaginary part negative) of the left half
		(fftsize/2-1).do({
			| i |
			re[i+(fftsize/2)] = re[(fftsize/2)-i];
			im[i+(fftsize/2)] = im[(fftsize/2)-i].neg;
		});

		// inverse fourier transformation: resulting imaginary part should be (very close to) all zeros
		re = ifft(re, im, tab);

		// re.real.normalize scales the result so it falls between 0 and 1.
		// Next, make sure to normalize the maximum volume to -3dB.
		result = re.real.normalize * ((-3).dbamp);

		// load the result in a buffer
		buffer = Buffer.loadCollection(s, result);

		// and return the buffer as result of the function
		buffer;
	};

	// calculate 2 wavetables per octave
	prepareAllBuffers = {
		| partials = #[ 1.01, 0.1722,  2.00, 0.0056,  2.99, 0.1609,  3.99, 0.0333,  5.00, 0.1157,
			5.99, 0.1149,  6.98, 0.0079, 7.98, 0.0620,  8.99, 0.0601,  9.99, 0.0104,
			10.98, 0.0134, 11.97, 0.0122, 12.99, 0.0058, 13.98, 0.0110, 14.98, 0.0029,
			15.97, 0.0045, 16.98, 0.0023, 17.98, 0.0010, 18.97, 0.0016, 19.96, 0.0021,
			20.96, 0.0008, 21.97, 0.0021, 22.96, 0.0001, 23.96, 0.0012, 24.95, 0.0003,
			25.97, 0.0002, 26.96, 0.0003, 27.95, 0.0002, 30.96, 0.0002, 32.94, 0.0002,
			34.96, 0.0001, 35.95, 0.0002, 37.93, 0.0001 ],
		  min_length=5,
		  spread = 0.1 |
		var buffers = [];
		var maxOctaves = 8;
		// prepare two wavetables per octave
		(0,1..8).do({
			| octave |
			var reference_note_1 = (octave*12);
			var reference_note_2 = (octave*12) + 6;
			buffers = buffers.add(prepareSingleBuffer.value(partials, min_length, spread, reference_note_1));
			buffers = buffers.add(prepareSingleBuffer.value(partials, min_length, spread, reference_note_2));
		});

		buffers;
	};

	// calculate a desired spectrum (note: envelopes and filters/resonators/fx are just as important in determining overall experience)
	xvals = (1,2..33).as(Array);
	spectrum = Signal.newClear(33).waveFill({
		| x, old, idx |
		var lookup;
		lookup = [ 14.2, 8.8, 7.3, 8, 5.7, 7, 6.8, 5.8, 8.7, 6.9, 3.2, 2.1, 4, 3, 1.8, 1.1, 2.5, 1.5]; // cello-esque if attack 0.3, no filter env, small reverb
		if ((idx < lookup.size), {lookup[idx]}, {0});

	}, start:1, end:0).as(Array);
	spectrum = [xvals, spectrum].lace;

	// prepare the wavetables using inverse FFT
	~buffers = prepareAllBuffers.value(spectrum,
		5 /* minimum length of buffer in seconds */,
		0.03 /* spread of partials during detuning */);

	s.sync;

	// start the fx synths
	reverb = Synth(\reverb, [
		\out, 0,
		\inbus, reverbbus,
		\mix, 0.1,
		\room, 0.5,
	],
	target:fxgroup2);

	chorus = Synth(\chorus, [
		\out, reverbbus,
		\inbus, chorusbus
	],
	target:fxgroup1);

	// create a composition
	// Pbind a new synth for each note.
	p = Pbind(
		\instrument, \padfilterenv,
		\out, reverbbus,
		\mynote, Pseq([Pseq((40,41..72), 1), Prand((48,49..72), 100)], inf),
		\myreference, (Pkey(\mynote)-(Pkey(\mynote)%6)),
		\referencefreq, Pkey(\myreference).midicps,
		\freq, Pkey(\mynote).midicps,
		\dur, Pseq([Pseq([0.5], 72-24), Prand((0.1,0.2..1.0), 100)], inf),
		\amp, 0.7,
		\buffer, Pfunc({ |ev| ~buffers[(ev[\myreference]/6).round(1)].bufnum; }),
		\attack, 0.3,
		\decay, 0.1,
		\sustain, 0.9,
		\release, 1.0,
		\filtergain, 0.3,
		\filtercutoff, 1000,
		\filterattack, 0.01,
		\filterdecay, 0.1,
		\filtersustain, 1.0,
		\filterrelease, 0.3,
		\filterresonance, 1.0,
		\glissando, 0,
		\vibratofreq, 3.0,
		\vibratodepth, 0.015,
		\group, sgroup);

	// Pmono creates only one synth, and updates its parameters. This allows e.g. for glissando's.
	q = Pmono(
		\padfilterenv,
		\out, reverbbus,
		\mynote, Pseq([40, 52], inf),
		\myreference, (Pkey(\mynote)-(Pkey(\mynote)%6)),
		\referencefreq, Pkey(\myreference).midicps,
		\freq, Pkey(\mynote).midicps,
		\dur, Pseq([2], inf),
		\amp, 0.7,
		\buffer, Pfunc({ |ev| ~buffers[(ev[\myreference]/6).round(1)].bufnum; }),
		\attack, 0.3,
		\decay, 0.1,
		\sustain, 0.9,
		\release, 1.0,
		\filtergain, 0.3,
		\filtercutoff, 1000,
		\filterattack, 0.01,
		\filterdecay, 0.1,
		\filtersustain, 1.0,
		\filterrelease, 0.3,
		\filterresonance, 1.0,
		\glissando, 0.5,
		\vibratofreq, 3.0,
		\vibratodepth, 0.015,
		\group, sgroup);

	c = Ppar([p,q], inf);
	c.play;

});
)

s.freeAll;
s.freeAllBuffers;
