<?php

/*
	NeuroElf remote web interface
	Server-side script (PHP) component

Version:  v0.9d
Build:    14071614
Date:     Jul-16 2014, 2:42 PM EST
Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
URL/Info: http://neuroelf.net/

Copyright (c) 2012 - 2014, Jochen Weber
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
	* Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	* Neither the name of Columbia University nor the
	  names of its contributors may be used to endorse or promote products
	  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* configuration */
$ne_folder = '/Applications/NeuroElf/NeuroElf_v10_5167/_files/cache';

/* no user cookie */
if (!array_key_exists('user', $_COOKIE) ||
	!array_key_exists('passwd', $_COOKIE) ||
	!array_key_exists('nesess', $_COOKIE)) {

	/* user and password given in URL-GET requests */
	if (array_key_exists('user', $_POST) &&
		array_key_exists('passwd', $_POST)) {

		/* authenticate */
		$tid = mt_rand(0, 16777215);
		$cid = mt_rand(0, 99999999);
		$tf = sprintf('%s/ner-%06x.nrt', $ne_folder, $tid);
		$fid = fopen($tf, 'a');
		fwrite($fid, $cid . '%auth%' . $_POST['user'] . '%' . $_POST['passwd']);
		fclose($fid);

		/* give it up to 15 seconds time */
		$c = 0;
		$of = sprintf('%s/ner-%06x-%s.nro', $ne_folder, $tid, $cid . '');
		$tc = '';
		while ($c < 1500) {

			/* see if file exists */
			if (file_exists($of)) {

				/* read in file and delete afterwards */
				$tc = file_get_contents($of);
				unlink($of);

				/* no more execution of loop */
				break;
			}

			/* pause */
			usleep(10000);
			$c++;
		}

		/* remove text file */
		if (file_exists($tf))
			unlink($tf);

		/* nothing read */
		if ($tc === '') {
			$tc = 'failed';
		}

		/* set cookie if OK */
		if (substr($tc, 0, 2) == 'ok') {
			$user = $_POST['user'];
			setcookie('user', $user);
			$passwd = $_POST['passwd'];
			setcookie('passwd', $passwd);
			setcookie('nesess', substr($tc, 2, 6));

		/* otherwise */
		} else {

			/* show login failed */
			readfile( __DIR__ . '/loginfailed.hti' );
			exit;
		}

	/* show login interface */
	} else {
		readfile( __DIR__ . '/login.hti' );
		exit;
	}

/* cookie present */
} else {
	$user = $_COOKIE['user'];
	$passwd = $_COOKIE['passwd'];
	$nesess = $_COOKIE['nesess'];
}

/* command requested */
if (array_key_exists('cmd', $_GET)) {

	/* get command */
	$command = $_GET['cmd'];

	/* write out new file */
	$tid = mt_rand(0, 16777215);
	$cid = mt_rand(0, 99999999);
	$tf = sprintf('%s/ner-%06x.nrt', $ne_folder, $tid);
	$fid = fopen($tf, 'a');
	fwrite($fid, $cid .'%'. $user .'%'. $passwd .'%'. $nesess .'%'. $command);
	fclose($fid);

	/* give it up to 15 seconds time */
	$c = 0;
	$of = sprintf('%s/ner-%06x-%s.nro', $ne_folder, $tid, $cid . '');
	$tc = '';
	while ($c < 1500) {

		/* see if file exists */
		if (file_exists($of)) {

			/* for images */
			if (preg_match('/^(corslice|sagslice|traslice|zoomslice|uishot)/i', $command) && (filesize($of) > 0)) {

				/* send image and exit */
				$tc = file_get_contents($of);
				unlink($of);

				/* try to detect filetype */
				if (substr($tc, 1, 3) == 'PNG')
					header('Content-type: image/png');
				else
					header('Content-type: image/jpeg');

			/* for text files */
			} else {

				/* read file into tc */
				$tc = file_get_contents($of);
				unlink($of);

				/* assume text header */
				header('Content-type: text/plain');
			}

			/* no more execution of loop */
			break;
		}

		/* pause */
		usleep(10000);
		$c++;
	}

	/* remove text file */
	if (file_exists($tf))
		unlink($tf);

	/* nothing read */
	if ((strlen($tc) == 0) ||
		(preg_match('/session_timeout/', $tc))) {

		/* for specific commands */
		if (preg_match('/^(corslice|sagslice|traslice)/i', $command)) {

			/* read default image */
			header('Content-type: image/png');
			readfile( __DIR__ . '/images/timeout.png');
			exit;

		/* otherwise */
		} else if (preg_match('/session_timeout/', $tc)) {
			setcookie('user', '', 1);
			setcookie('passwd', '', 1);
			setcookie('nesess', '', 1);
		}
	}

	print $tc;

} else {

	/* another but the main form requested? */
	if (array_key_exists('form', $_GET) &&
		preg_match('/^[a-z_]+$/i', $_GET['form']) &&
		file_exists( __DIR__ . '/' . $_GET['form'] . '.hti' ))
		$form = $_GET['form'];
	else
		$form = 'ne_remote';

	/* show interface/form */
	readfile( __DIR__ . '/' . $form . '.hti' );
}
?>
