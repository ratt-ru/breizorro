from bokeh.io import curdoc
from bokeh.layouts import row
from bokeh.models import (
    CustomJS,
    DataTable,
    HoverTool,
    LabelSet,
    TableColumn,
)
from bokeh.plotting import ColumnDataSource, figure, output_file, show

from breizorro.utils import fitsInfo, format_source_coordinates

curdoc().theme = "caliber"


def display(imagename, mask_image, outcatalog, source_list):
    """
    Display the image with overlaid source positions and associated catalog information.

    Parameters
    ----------
    imagename : str
        The name or path of the FITS file to be displayed.

    outcatalog : bool
        If True, overlay the detected sources on the image using the provided source list.

    source_list : list of tuples
        A list containing information about detected sources, where each entry is a tuple
        with the following format:
        (RA, 'RA DEC RA_err DEC_err I I_err I_peak I_peak_error Emaj_s Emin_s PA_d', flag).

        - RA : float
            Right Ascension in degrees.
        - DEC : float
            Declination in degrees.
        - RA_err : float
            Right Ascension uncertainty in degrees.
        - DEC_err : float
            Declination uncertainty in degrees.
        - I : float
            Intensity or flux measurement.
        - I_err : float
            Error in the integrated intensity measurement.
        - I_peak : float
            Peak intensity.
        - I_peak_error : float
            Error in the peak intensity measurement.
        - Emaj_s : float
            Major axis error (in arcseconds).
        - Emin_s : float
            Minor axis error (in arcseconds).
        - PA_d : float
            Position angle (in degrees).
        - flag : int
            A flag indicating some property of the source (e.g., detection confidence).

    Returns
    -------
    None
        Displays/Saves an image with overlaid catalog scatter points and mask image.
    """
    # Get fits information
    fitsinfo = fitsInfo(imagename)
    # Origin coordinates
    origin_ra, origin_dec = fitsinfo["centre"]
    # Ensure RA is always positive (in degrees, 0 to 360)
    origin_ra = origin_ra % 360
    if outcatalog:
        # Extracting data from source_list
        # x_coords = [float(d[1].split(" ")[0]) for d in source_list]
        y_coords = [float(d[1].split(" ")[1]) for d in source_list]
        labels = [
            f"{format_source_coordinates(float(d[1].split(' ')[0]), float(d[1].split(' ')[1]))}" for d in source_list
        ]

        # Assuming `source_list` is already populated with your data
        # Example source_list: [(ra, 'ra dec ra_err dec_err i i_err i_peak i_peak_error emaj_s emin_s pa_d', flag)]
        # Parse the source_list and split each string into its components
        data = {
            "name": [f"src{i}" for i in range(len(source_list))],
            "ra_deg": [float(d[1].split(" ")[0]) for d in source_list],
            "dec_deg": [float(d[1].split(" ")[1]) for d in source_list],
            "ra_err": [float(d[1].split(" ")[2]) for d in source_list],
            "dec_err": [float(d[1].split(" ")[3]) for d in source_list],
            "i": [float(d[1].split(" ")[4]) for d in source_list],
            "error": [float(d[1].split(" ")[5]) for d in source_list],
            "i_peak": [float(d[1].split(" ")[6]) for d in source_list],
            "i_peak_error": [float(d[1].split(" ")[7]) for d in source_list],
            "emaj_s": [float(d[1].split(" ")[8]) for d in source_list],
            "emin_s": [float(d[1].split(" ")[9]) for d in source_list],
            "pa_d": [float(d[1].split(" ")[10]) for d in source_list],
        }

        # Format RA and DEC to hh:mm:ss and dd:mm:ss
        formatted_coords = [format_source_coordinates(ra, dec) for ra, dec in zip(data["ra_deg"], data["dec_deg"])]
        formatted_RA = [coord[0] for coord in formatted_coords]
        formatted_DEC = [coord[1] for coord in formatted_coords]

        # Add formatted coordinates to the data source
        data["RA"] = formatted_RA
        data["DEC"] = formatted_DEC

        # Create a ColumnDataSource for the table
        table_source = ColumnDataSource(data=data)
        from bokeh.models import NumberFormatter

        # Define the number formatter for decimal and scientific notation
        decimal_formatter = NumberFormatter(format="0.000")
        scientific_formatter = NumberFormatter(format="0.000e")
        # Define the columns for the DataTable
        columns = [
            TableColumn(field="name", title="name"),
            TableColumn(field="ra_deg", title="ra (deg)", formatter=decimal_formatter),
            TableColumn(field="dec_deg", title="dec (deg)", formatter=decimal_formatter),
            TableColumn(field="ra_err", title="ra_err (deg)", formatter=scientific_formatter),
            TableColumn(field="dec_err", title="dec_err (deg)", formatter=scientific_formatter),
            TableColumn(field="i", title="i (Jy)", formatter=scientific_formatter),
            TableColumn(field="error", title="i_err (Jy)", formatter=scientific_formatter),
            TableColumn(field="i_peak", title="i_peak (Jy)", formatter=scientific_formatter),
            TableColumn(field="i_peak_error", title="i_peak_error (Jy)", formatter=scientific_formatter),
            TableColumn(field="emaj_s", title="emaj_s (arcsec)"),
            TableColumn(field="emin_s", title="emin_s (arcsec)"),
            TableColumn(field="pa_d", title="pa_d (deg)"),
        ]

        # Create the DataTable widget
        data_table = DataTable(source=table_source, columns=columns, width=600, height=600)

        # Create the plot
        p = figure(
            title="Breizorro Source Catalog",
            x_axis_label="Right Ascension",
            y_axis_label="Declination",
            # x_range=(max(x_coords), min(x_coords)),
            y_range=(min(y_coords), max(y_coords)),
            match_aspect=True,
            tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
        )
        scatter_renderer = p.scatter(
            "ra_deg",
            "dec_deg",
            size=6,
            color="red",
            source=table_source,
            legend_label="Detection",
            level="overlay",
        )
        # Add labels to the scatter points (optional, can hide later as needed)
        labels = LabelSet(
            x="ra_deg",
            y="dec_deg",
            text="name",
            source=table_source,
            text_font_size="10pt",
            text_baseline="middle",
            text_align="center",
        )

        # Add hover tool for scatter points
        hover = HoverTool()
        hover.tooltips = [
            ("Name", "@name"),
            ("RA", "@RA"),
            ("DEC", "@DEC"),
            ("Flux", "@i"),
            ("Peak", "@i_peak"),
        ]
        hover.renderers = [scatter_renderer]
        p.add_tools(hover)
        # Enable legend click to hide/show scatter points
        p.legend.click_policy = "hide"
        p.legend.label_text_color = "white"  # Set the legend text color to white
        p.x_range.flipped = True
        p.title.align = "center"
        p.add_layout(labels)
        p.grid.grid_line_width = 0.5
        # from bokeh.models import DatetimeTickFormatter
        # p.xaxis.formatter = DatetimeTickFormatter(hours=["%Hh%Mm%Ss"])
        # p.yaxis.formatter = DatetimeTickFormatter(minutes=["%Dd%Mm%Ss"])
        from bokeh.models import CustomJSTickFormatter

        # Formatter for RA axis (converting degrees to hh:mm:ss)
        p.xaxis.formatter = CustomJSTickFormatter(
            code="""
            // Convert RA from degrees to hours, minutes, and seconds
            let hours = Math.floor(tick / 15);  // 15 degrees per hour
            let minutes = Math.floor((tick % 15) * 4);
            let seconds = Math.floor((((tick % 15) * 4) - minutes) * 60);
            return hours + "h " + minutes + "m " + seconds + "s";
        """
        )

        # Formatter for Dec axis (converting degrees to dd:mm:ss)
        p.yaxis.formatter = CustomJSTickFormatter(
            code="""
            // Convert Dec from degrees to degrees, minutes, and seconds
            let degrees = Math.floor(Math.abs(tick));
            let minutes = Math.floor((Math.abs(tick) - degrees) * 60);
            let seconds = Math.floor((((Math.abs(tick) - degrees) * 60) - minutes) * 60);
            let sign = tick < 0 ? '-' : '';
            return sign + degrees + "° " + minutes + "' " + seconds + '"';
        """
        )
        # JavaScript callback to select the corresponding row in the table when hovering over a point
        callback = CustomJS(
            args=dict(source=table_source),
            code="""
            var indices = cb_data.index['1d'].indices;  // Get the index of the hovered point
            source.selected.indices = indices;  // Select the corresponding row in the table
        """,
        )
        # Connect the hover tool to the callback
        hover.callback = callback
        layout = row(p, data_table, sizing_mode="stretch_both")
        output_file("breizorro.html", title="Mask Editor")
        # export_png(p, filename="breizorro.png")
        show(layout)
