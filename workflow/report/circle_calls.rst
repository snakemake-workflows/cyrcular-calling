{% set scenario_events = snakemake.config['filter']['fdr-control']['events'][snakemake.wildcards.event]['varlociraptor'] %}
{% set fdr = snakemake.config['filter']['fdr-control']['threshold'] %}
{% set mode = snakemake.config['filter']['fdr-control']['mode'].split(-) %}
Circle calls for the **event {{snakemake.wildcards.event}}**, as defined in ``config/config.yaml`` (encompassing the varlociraptor scenario events: {{ scenario_events|join(", ") }} ).

{% if "global" in mode %}
Calling was done in `**global** false discovery rate <https://en.wikipedia.org/wiki/False_discovery_rate>`_ mode, meaning that the global expected fraction of false positives has to stay below the configured threshold of {{ fdr }}.
{% else %}
Calling was done in `**local** false discovery rate <https://en.wikipedia.org/wiki/False_discovery_rate>`_ mode, meaning that the false discovery probability for each individual variant has to stay below the configured threshold of {{ fdr }}.
{% endif %}
{% if "smart" in mode %}
As this was done in **smart mode**, the probability of a false discovery is calculated as the sum of the probabilities that the circle is absent or an artifact and a circle is only selected, if the sum of the probabilities of the varlociraptor scenario events {{ scenario_events|join(", ") }} is greater than any of the other event probabilities.
{% else %}
As this was done in **strict mode**, the probablity of a false discovery is calculated as 1 minus the sum of the varlociraptor scenario events {{ scenario_events|join(", ")}}.
{% endif %}
For more details, see the `varlociraptor documentation on filtering <https://varlociraptor.github.io/docs/filtering/>`_.