# Generating a flycode FASTA database for LC-MS/MS
```SQL
-- Grouping by barcode, then variants and concatenating corresponding flycodes
CREATE OR REPLACE VIEW variant_flycodes AS

-- flycodes
WITH all_variants AS (
  SELECT
    barcode_id,
    barcode as flycode,
    "position" -1 AS "position",
    variant_aa,
    reference_aa || "position" -1 || variant_aa AS aa_change
  FROM variants
  WHERE
    variant_type = 'change' AND
    "position" -1 in (5, 9, 12, 47, 48, 50, 54, 67, 117, 120, 127, 158, 159, 166, 169, 193) AND
    variant_aa != '*'
  GROUP BY
    barcode_id, "position", barcode, reference_aa, variant_aa
),

-- only one amino acid change per barcode (flycode)
one_aa_change AS (
  SELECT barcode_id,
  FROM variants
  WHERE variant_type = 'change'
  GROUP BY barcode_id
  HAVING count("position") = 1
)

SELECT
  aa_change,
  "position",
  variant_aa,
  string_agg(flycode, '') AS flycodes
FROM all_variants
JOIN one_aa_change using(barcode_id)
GROUP BY
  "position", variant_aa, aa_change
ORDER BY
  "position" ASC, variant_aa ASC;

-- Writing FASTA file
COPY (
  SELECT '>KDELR=' || aa_change || chr(10) || flycodes,
  FROM variant_flycodes
) TO 'flycode_db.fasta' (HEADER FALSE, DELIMITER '', QUOTE '');
```