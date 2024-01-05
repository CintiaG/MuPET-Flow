## Help

Sorry, we are in construction...

### Known issues

* App works deletes `.fcs` and `.FCS` extensions from files names, and removes everything after a space to generate the sample name. If after this previous process duplicated names are generated, the app will crash, thus it is better to name the files with unique names
* Sometimes it crashes when uploading different samples, so it is better to close and re-run the app.
* Notice that when modifying the smoothing and window, it recalculates and deselect already chosen peaks.
* When files are large and exceed the maximum size, prior running the app you can set in your console: `options(shiny.maxRequestSize = 10 * 1024^2)`.
* It will give an error if window (span) is too small.

### Potential improvements
* Add safe code to avoid downloading empty documents.
* Add code to download individual histograms and linear regression plots.
* Add code to download regression summary as text.
* Add code to exclude samples from final summary histograms.
* Add code to download all detected peaks when more than two are detected.
* Generate new version of script to resume analysis using previously determined smoothing and window.
* Explain in help what is the plotted r and the pvalue.
