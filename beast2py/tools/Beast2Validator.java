import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.pkgmgmt.PackageManager;

import java.io.File;

/**
 * Headless BEAST2 XML validator.
 *
 * Uses beast.base.parser.XMLParser to parse and validate an XML file
 * without requiring JavaFX (unlike beastfx.app.beast.BeastMain).
 *
 * Exit codes:
 *   0 = XML is valid
 *   1 = XML is invalid (parse error, missing class, etc.)
 *   2 = Usage error
 */
public class Beast2Validator {

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: Beast2Validator <xml-file> [<xml-file> ...]");
            System.exit(2);
        }

        // Load all installed BEAST2 packages
        try {
            PackageManager.loadExternalJars();
        } catch (Exception e) {
            System.err.println("WARNING: Failed to load BEAST2 packages: " + e.getMessage());
            // Continue anyway - classes may be on the classpath directly
        }

        boolean allValid = true;
        for (String xmlPath : args) {
            File xmlFile = new File(xmlPath);
            if (!xmlFile.exists()) {
                System.err.println("INVALID: File not found: " + xmlPath);
                allValid = false;
                continue;
            }

            try {
                XMLParser parser = new XMLParser();
                parser.parseFile(xmlFile);
                System.out.println("VALID: " + xmlPath);
            } catch (XMLParserException e) {
                System.err.println("INVALID: " + xmlPath);
                System.err.println("  " + e.getMessage());
                allValid = false;
            } catch (Exception e) {
                System.err.println("INVALID: " + xmlPath);
                System.err.println("  " + e.getClass().getSimpleName() + ": " + e.getMessage());
                if (e.getCause() != null) {
                    System.err.println("  Caused by: " + e.getCause().getClass().getSimpleName()
                        + ": " + e.getCause().getMessage());
                }
                allValid = false;
            }
        }

        System.exit(allValid ? 0 : 1);
    }
}
