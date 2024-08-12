from pathlib import Path
from extractor.article import Article
from extractor.parsehtml import Parser

article_AAA = Article(title = "Automated grading systems for programming assignments: A literature review", authors = ("H Aldriye", "A Alkhalaf", "M Alkhalaf"), authors_partial=False, year=2019, abstract="… open-source and online automated grading system for students Java and C++ programming assignment… Web-CAT is a state of the art in automated grading system, and many instructors …", link=None, PDF="https://pdfs.semanticscholar.org/3ec5/6c7bc0338717acca9d5e519c0a27ca5ce70d.pdf")
article_IAK = Article(title="Review of recent systems for automatic assessment of programming assignments", authors = ("P Ihantola", "T Ahoniemi", "V Karavirta"), authors_partial=True, year=2010, abstract="… Assessment provides both means to guide student learning and … for Java or have support for Java. This fits well with the trend of … system developers to make their systems open source …", link="https://dl.acm.org/doi/abs/10.1145/1930464.1930480", PDF="https://www.researchgate.net/profile/Petri-Ihantola/publication/216714976_Review_of_recent_systems_for_automatic_assessment_of_programming_assignments/links/0deec5182391eb2954000000/Review-of-recent-systems-for-automatic-assessment-of-programming-assignments.pdf")
article_LKS23 = Article(title="Pedagogic Exploration Into Adapting Automated Writing Evaluation and Peer Review Integrated Feedback Into Large-Sized University Writing Classes", authors = ("WY Li", "K Kau", "YJ Shiung"), authors_partial=False, year=2023, abstract="… Grammarly for feedback and reflect upon their revisions for each assignment, and (3) establishing a Grammarly-peer review integrated feedback mode (IFM) during post-writing peer …", link="https://journals.sagepub.com/doi/abs/10.1177/21582440231209087", PDF="https://journals.sagepub.com/doi/pdf/10.1177/21582440231209087")
article_S21 = Article(title="A survey and evaluation of browser fingerprinting techniques", authors = ("W Segers", ), authors_partial=False, year=2021, abstract="… ) or to conduct automated tests on a website in development [29]. … Each fingerprinting vector in the fingerprint library adheres to … Where Flash and Java could probe to fonts and webcams …", link=None, PDF="https://documentserver.uhasselt.be/bitstream/1942/35316/1/fae24997-cffe-46e1-a2b6-2d8a542c18be.pdf")
article_SW13 = Article(title="Enabling effective synoptic assessment via algorithmic constitution of review panels", authors = ("ND de Silva", "SM Weerawarana"), authors_partial=True, year=2013, abstract="… Java is taught as an object oriented programming language … streamed into their respective open source codebases. Some … of expert evaluation panels by automating the panel …", link="https://ieeexplore.ieee.org/abstract/document/6654543/", PDF="https://www.researchgate.net/profile/Nisansa-De-Silva/publication/261164180_Enabling_effective_synoptic_assessment_via_algorithmic_constitution_of_review_panels/links/5948ffb6458515db1fdb3768/Enabling-effective-synoptic-assessment-via-algorithmic-constitution-of-review-panels.pdf")
article_KHR_19 = Article(title="A Survey on Automated Answer Paper Evaluation and Text Extraction Techniques", authors = ("N Kalaiselvi", "M Hariharan", "S Ramesh", "A Vignesh"), authors_partial=False, year=2019, abstract="", link=None)
article_DMP21 = Article(title="Automatic question generation and answer assessment: a survey", authors = ("B Das", "M Majumder", "S Phadikar"), authors_partial=True, year=2021, abstract="… facilitates learners to learn anything, anytime… survey of automatic question generation and assessment strategies from textual and pictorial learning resources. The purpose of this survey …", link="https://telrp.springeropen.com/articles/10.1186/s41039-021-00151-1")
article_B11 = Article(title = "Automatic Assessment of Java Programming Patterns for Novices", authors = ("A Bagini", ), authors_partial=False, year = 2011, abstract = "", link = None)
article_SE12 = Article(title = "An intelligent assessment tool for students' Java submissions in introductory programming courses", authors = ("FA Shamsi", "A Elnagar"), authors_partial=False, year = 2012, abstract = "", link = "https://www.scirp.org/html/6-9601091_17557.htm")
article_SGK = Article(title = "The vital role of community in open source software development: A framework for assessment and ranking", authors = ("J Singh", "A Gupta", "P Kanwal"), authors_partial=False, year = None, abstract = "… teams, identifies core FOSS evaluation criteria and introduces a fully automated method for … primary data for evaluation. To demonstrate the effectiveness of Open Source Quality Radar, …", link = "https://onlinelibrary.wiley.com/doi/abs/10.1002/smr.2643")

def test_entries():
    Parser.article_entries("")
    with open("tests/resources/res1.html", "r") as f:
        html = f.read()
    entries = Parser.article_entries(html)
    assert len(entries) == 10

def test_article_AAA():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="-uN63_al6jUJ" data-did="-uN63_al6jUJ" data-lid="IWHjjKOFINEC" data-aid="-uN63_al6jUJ" data-rp="2">
        <div class="gs_ggs gs_fl">
            <div class="gs_ggsd">
                <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://pdfs.semanticscholar.org/3ec5/6c7bc0338717acca9d5e519c0a27ca5ce70d.pdf" data-clk="hl=en&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=2&amp;d=3885100108290384890&amp;ei=ODyEZZSGC5iEy9YPgLGi6AI" data-clk-atid="-uN63_al6jUJ"><span class="gs_ctg2">[PDF]</span> semanticscholar.org</a></div>
            </div>
        </div>
        <div class="gs_ri">
            <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><span class="gs_ctc"><span class="gs_ct1">[PDF]</span><span class="gs_ct2">[PDF]</span></span> <a id="-uN63_al6jUJ" href="https://pdfs.semanticscholar.org/3ec5/6c7bc0338717acca9d5e519c0a27ca5ce70d.pdf" data-clk="hl=en&amp;sa=T&amp;oi=ggp&amp;ct=res&amp;cd=2&amp;d=3885100108290384890&amp;ei=ODyEZZSGC5iEy9YPgLGi6AI" data-clk-atid="-uN63_al6jUJ"><b>Automated grading </b>systems for programming assignments: A literature <b>review</b></a></h3>
            <div class="gs_a">H Aldriye, A Alkhalaf, M Alkhalaf&nbsp;- International Journal of&nbsp;…, 2019 - pdfs.semanticscholar.org</div>
            <div class="gs_rs">… <b>open</b>-<b>source</b> and online <b>automated</b> <b>grading</b> system for students <b>Java</b> and C++ programming <br>
                <b>assignment</b>… Web-CAT is a state of the art in <b>automated</b> <b>grading</b> system, and many instructors …
            </div>
            <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                </svg><span>Cite</span></a> <a href="https://scholar.google.com/scholar?cites=3885100108290384890&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=en">Cited by 38</a> <a href="https://scholar.google.com/scholar?q=related:-uN63_al6jUJ:scholar.google.com/&amp;scioq=(intitle:survey+OR+intitle:review)++AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:CBA)++AND+(intext:java)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)++AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+&amp;hl=en&amp;as_sdt=0,5">Related articles</a> <a href="https://scholar.google.com/scholar?cluster=3885100108290384890&amp;hl=en&amp;as_sdt=0,5" class="gs_nph">All 2 versions</a> <a href="javascript:void(0)" title="More" class="gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                </svg></a> <a href="https://scholar.googleusercontent.com/scholar?q=cache:-uN63_al6jUJ:scholar.google.com/+(intitle:survey+OR+intitle:review)++AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:CBA)++AND+(intext:java)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)++AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+&amp;hl=en&amp;as_sdt=0,5" class="gs_or_nvi">View as HTML</a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                </svg></a></div>
        </div>
    </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_AAA

def test_article_IAK():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="Yf7ALe8h_DEJ" data-did="Yf7ALe8h_DEJ" data-lid="" data-aid="Yf7ALe8h_DEJ" data-rp="3">
        <div class="gs_ggs gs_fl">
            <div class="gs_ggsd">
                <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://www.researchgate.net/profile/Petri-Ihantola/publication/216714976_Review_of_recent_systems_for_automatic_assessment_of_programming_assignments/links/0deec5182391eb2954000000/Review-of-recent-systems-for-automatic-assessment-of-programming-assignments.pdf" data-clk="hl=en&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=3&amp;d=3601791113138077281&amp;ei=ODyEZZSGC5iEy9YPgLGi6AI" data-clk-atid="Yf7ALe8h_DEJ"><span class="gs_ctg2">[PDF]</span> researchgate.net</a></div>
            </div>
        </div>
        <div class="gs_ri">
            <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><a id="Yf7ALe8h_DEJ" href="https://dl.acm.org/doi/abs/10.1145/1930464.1930480" data-clk="hl=en&amp;sa=T&amp;ct=res&amp;cd=3&amp;d=3601791113138077281&amp;ei=ODyEZZSGC5iEy9YPgLGi6AI" data-clk-atid="Yf7ALe8h_DEJ"><b>Review </b>of recent systems for <b>automatic assessment </b>of programming assignments</a></h3>
            <div class="gs_a"><a href="https://scholar.google.com/citations?user=8rxsMsYAAAAJ&amp;hl=en&amp;oi=sra">P Ihantola</a>, T Ahoniemi, <a href="https://scholar.google.com/citations?user=867GoxcAAAAJ&amp;hl=en&amp;oi=sra">V Karavirta</a>…&nbsp;- Proceedings of the 10th&nbsp;…, 2010 - dl.acm.org</div>
            <div class="gs_rs">… <b>Assessment</b> provides both means to guide <b>student</b> learning and … for <b>Java</b> or have support <br>
                for <b>Java</b>. This fits well with the trend of … system developers to make their systems <b>open</b> <b>source</b> …</div>
            <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                </svg><span>Cite</span></a> <a href="https://scholar.google.com/scholar?cites=3601791113138077281&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=en">Cited by 637</a> <a href="https://scholar.google.com/scholar?q=related:Yf7ALe8h_DEJ:scholar.google.com/&amp;scioq=(intitle:survey+OR+intitle:review)++AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:CBA)++AND+(intext:java)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)++AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+&amp;hl=en&amp;as_sdt=0,5">Related articles</a> <a href="https://scholar.google.com/scholar?cluster=3601791113138077281&amp;hl=en&amp;as_sdt=0,5" class="gs_nph">All 7 versions</a> <a href="javascript:void(0)" title="More" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                </svg></a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                    <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                </svg></a></div>
        </div>
    </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_IAK

def test_article_S21():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="6B-P5fW1z6oJ" data-did="6B-P5fW1z6oJ" data-lid="" data-aid="6B-P5fW1z6oJ" data-rp="144">
              <div class="gs_ggs gs_fl">
                <div class="gs_ggsd">
                  <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://documentserver.uhasselt.be/bitstream/1942/35316/1/fae24997-cffe-46e1-a2b6-2d8a542c18be.pdf" data-clk="hl=en&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=144&amp;d=12308256374349832168&amp;ei=jsOEZeLFMIOGy9YPsbew6A0" data-clk-atid="6B-P5fW1z6oJ"><span class="gs_ctg2">[PDF]</span> uhasselt.be</a></div>
                </div>
              </div>
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><span class="gs_ctc"><span class="gs_ct1">[PDF]</span><span class="gs_ct2">[PDF]</span></span> <a id="6B-P5fW1z6oJ" href="https://documentserver.uhasselt.be/bitstream/1942/35316/1/fae24997-cffe-46e1-a2b6-2d8a542c18be.pdf" data-clk="hl=en&amp;sa=T&amp;oi=ggp&amp;ct=res&amp;cd=144&amp;d=12308256374349832168&amp;ei=jsOEZeLFMIOGy9YPsbew6A0" data-clk-atid="6B-P5fW1z6oJ">A <b>survey </b>and <b>evaluation </b>of browser fingerprinting techniques</a></h3>
                <div class="gs_a">W Segers - 2021 - documentserver.uhasselt.be</div>
                <div class="gs_rs">… ) or to conduct <b>automated</b> tests on a website in development [29]. … Each fingerprinting vector <br>
                  in the fingerprint <b>library</b> adheres to … Where Flash and <b>Java</b> could probe to fonts and webcams …</div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Cite</span></a> <a href="https://scholar.google.com/scholar?q=related:6B-P5fW1z6oJ:scholar.google.com/&amp;scioq=(intitle:survey+OR+intitle:review)++AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:CBA)++AND+(intext:java)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)++AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+&amp;hl=en&amp;as_sdt=0,5">Related articles</a> <a href="javascript:void(0)" title="More" class="gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="https://scholar.googleusercontent.com/scholar?q=cache:6B-P5fW1z6oJ:scholar.google.com/+(intitle:survey+OR+intitle:review)++AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:CBA)++AND+(intext:java)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)++AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+&amp;hl=en&amp;as_sdt=0,5" class="gs_or_nvi">View as HTML</a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_S21

def test_article_SW13():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="6CyDz1cR6nQJ" data-did="6CyDz1cR6nQJ" data-lid="" data-aid="6CyDz1cR6nQJ" data-rp="97">
              <div class="gs_ggs gs_fl">
                <div class="gs_ggsd">
                  <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://www.researchgate.net/profile/Nisansa-De-Silva/publication/261164180_Enabling_effective_synoptic_assessment_via_algorithmic_constitution_of_review_panels/links/5948ffb6458515db1fdb3768/Enabling-effective-synoptic-assessment-via-algorithmic-constitution-of-review-panels.pdf" data-clk="hl=en&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=97&amp;d=8424565121791241448&amp;ei=gcOEZaTLBpLGy9YPl_WhiAY" data-clk-atid="6CyDz1cR6nQJ"><span class="gs_ctg2">[PDF]</span> researchgate.net</a></div>
                </div>
              </div>
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><a id="6CyDz1cR6nQJ" href="https://ieeexplore.ieee.org/abstract/document/6654543/" data-clk="hl=en&amp;sa=T&amp;ct=res&amp;cd=97&amp;d=8424565121791241448&amp;ei=gcOEZaTLBpLGy9YPl_WhiAY" data-clk-atid="6CyDz1cR6nQJ">Enabling effective synoptic <b>assessment </b>via algorithmic constitution of <b>review </b>panels</a></h3>
                <div class="gs_a"><a href="https://scholar.google.com/citations?user=OgPD66oAAAAJ&amp;hl=en&amp;oi=sra">ND de Silva</a>, SM Weerawarana…&nbsp;- …&nbsp;Teaching, <b>Assessment</b>&nbsp;…, 2013 - ieeexplore.ieee.org</div>
                <div class="gs_rs">… <b>Java</b> is taught as an object oriented programming language … streamed into their respective <br>
                  <b>open</b> <b>source</b> codebases. Some … of expert <b>evaluation</b> panels by <b>automating</b> the panel …
                </div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Cite</span></a> <a href="https://scholar.google.com/scholar?cites=8424565121791241448&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=en">Cited by 2</a> <a href="https://scholar.google.com/scholar?q=related:6CyDz1cR6nQJ:scholar.google.com/&amp;scioq=(intitle:survey+OR+intitle:review)++AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:CBA)++AND+(intext:java)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)++AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+&amp;hl=en&amp;as_sdt=0,5">Related articles</a> <a href="https://scholar.google.com/scholar?cluster=8424565121791241448&amp;hl=en&amp;as_sdt=0,5" class="gs_nph">All 7 versions</a> <a href="javascript:void(0)" title="More" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_SW13

def test_article_LKS23():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="avP9Q7KopP0J" data-did="avP9Q7KopP0J" data-lid="" data-aid="avP9Q7KopP0J" data-rp="155">
              <div class="gs_ggs gs_fl">
                <div class="gs_ggsd">
                  <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://journals.sagepub.com/doi/pdf/10.1177/21582440231209087" data-clk="hl=en&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=155&amp;d=18276918671374676842&amp;ei=jsOEZeLFMIOGy9YPsbew6A0" data-clk-atid="avP9Q7KopP0J"><span class="gs_ctg2">[PDF]</span> sagepub.com</a><a href="https://scholar.google.com/scholar?output=instlink&amp;q=info:avP9Q7KopP0J:scholar.google.com/&amp;hl=en&amp;as_sdt=0,5&amp;scillfp=4701059090905929268&amp;oi=lle">Full View</a></div>
                </div>
              </div>
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><a id="avP9Q7KopP0J" href="https://journals.sagepub.com/doi/abs/10.1177/21582440231209087" data-clk="hl=en&amp;sa=T&amp;ct=res&amp;cd=155&amp;d=18276918671374676842&amp;ei=jsOEZeLFMIOGy9YPsbew6A0" data-clk-atid="avP9Q7KopP0J">Pedagogic Exploration Into Adapting <b>Automated </b>Writing <b>Evaluation </b>and Peer <b>Review </b>Integrated Feedback Into Large-Sized University Writing Classes</a></h3>
                <div class="gs_a">WY Li, K Kau, YJ Shiung&nbsp;- SAGE Open, 2023 - journals.sagepub.com</div>
                <div class="gs_rs">… Grammarly for feedback and reflect upon their revisions for each <b>assignment</b>, and (3) <br>
                  establishing a Grammarly-peer <b>review</b> integrated feedback mode (IFM) during post-writing peer …</div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Cite</span></a> <a href="javascript:void(0)" title="More" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_LKS23

def test_article_KHR19():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="_HiGHDN2BzgJ" data-did="_HiGHDN2BzgJ" data-lid="" data-aid="_HiGHDN2BzgJ" data-rp="46">
    <div class="gs_ri">
      <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><span class="gs_ctu"><span class="gs_ct1">[CITATION]</span><span class="gs_ct2">[C]</span></span> <span id="_HiGHDN2BzgJ">A <b>Survey </b>on <b>Automated </b>Answer Paper <b>Evaluation </b>and Text Extraction Techniques</span></h3>
      <div class="gs_a">N Kalaiselvi, M Hariharan, <a href="https://scholar.google.com/citations?user=leik1dQAAAAJ&amp;hl=fr&amp;num=20&amp;oi=sra">S Ramesh</a>, A Vignesh - 2019</div>
      <div class="gs_rs"></div>
      <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
            <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
          </svg><span class="gs_or_btn_lbl">Enregistrer</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
            <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
          </svg><span>Citer</span></a> <a href="https://scholar.google.com/scholar?q=related:_HiGHDN2BzgJ:scholar.google.com/&amp;scioq=(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment)+(intitle:survey+OR+intitle:review)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)+AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)&amp;hl=fr&amp;num=20&amp;as_sdt=0,5">Autres articles</a> <a href="javascript:void(0)" title="Plus" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
            <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
          </svg></a> <a href="javascript:void(0)" title="Moins" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
            <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
          </svg></a></div>
    </div>
  </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_KHR_19

def test_article_DMP21():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="4oe-qecal9sJ" data-did="4oe-qecal9sJ" data-lid="" data-aid="4oe-qecal9sJ" data-rp="28">
              <div class="gs_ggs gs_fl">
                <div class="gs_ggsd">
                  <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://telrp.springeropen.com/articles/10.1186/s41039-021-00151-1" data-clk="hl=fr&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=28&amp;d=15823145398077130722&amp;ei=uF6MZfzJGJesy9YPv4Ok2Ao" data-clk-atid="4oe-qecal9sJ"><span class="gs_ctg2">[HTML]</span> springeropen.com</a><a href="https://scholar.google.com/scholar?output=instlink&amp;q=info:4oe-qecal9sJ:scholar.google.com/&amp;hl=fr&amp;num=20&amp;as_sdt=0,5&amp;scillfp=11541094906280355642&amp;oi=lle">Full View</a></div>
                </div>
              </div>
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><span class="gs_ctc"><span class="gs_ct1">[HTML]</span><span class="gs_ct2">[HTML]</span></span> <a id="4oe-qecal9sJ" href="https://telrp.springeropen.com/articles/10.1186/s41039-021-00151-1" data-clk="hl=fr&amp;sa=T&amp;oi=ggp&amp;ct=res&amp;cd=28&amp;d=15823145398077130722&amp;ei=uF6MZfzJGJesy9YPv4Ok2Ao" data-clk-atid="4oe-qecal9sJ"><b>Automatic </b>question generation and answer <b>assessment</b>: a <b>survey</b></a></h3>
                <div class="gs_a"><a href="https://scholar.google.com/citations?user=4_MPSs4AAAAJ&amp;hl=fr&amp;num=20&amp;oi=sra">B Das</a>, <a href="https://scholar.google.com/citations?user=WCp7DmIAAAAJ&amp;hl=fr&amp;num=20&amp;oi=sra">M Majumder</a>, <a href="https://scholar.google.com/citations?user=WIFP2yIAAAAJ&amp;hl=fr&amp;num=20&amp;oi=sra">S Phadikar</a>…&nbsp;- Research and Practice&nbsp;…, 2021 - telrp.springeropen.com</div>
                <div class="gs_rs">… facilitates learners to <b>learn</b> anything, anytime… <b>survey</b> of <b>automatic</b> question generation and <br>
                  <b>assessment</b> strategies from textual and pictorial learning resources. The purpose of this <b>survey</b> …
                </div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Enregistrer</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Citer</span></a> <a href="https://scholar.google.com/scholar?cites=15823145398077130722&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=fr&amp;num=20">Cité 48&nbsp;fois</a> <a href="https://scholar.google.com/scholar?q=related:4oe-qecal9sJ:scholar.google.com/&amp;scioq=(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment)+(intitle:survey+OR+intitle:review)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)+AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)&amp;hl=fr&amp;num=20&amp;as_sdt=0,5">Autres articles</a> <a href="https://scholar.google.com/scholar?cluster=15823145398077130722&amp;hl=fr&amp;num=20&amp;as_sdt=0,5" class="gs_nph">Les 9 versions</a> <a href="javascript:void(0)" title="Plus" class="gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="https://scholar.google.com/scholar?output=instlink&amp;q=info:4oe-qecal9sJ:scholar.google.com/&amp;hl=fr&amp;num=20&amp;as_sdt=0,5&amp;scillfp=4252293008347341383&amp;oi=llo" class="gs_or_nvi">Full View</a> <a href="https://scholar.googleusercontent.com/scholar?q=cache:4oe-qecal9sJ:scholar.google.com/+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment)+(intitle:survey+OR+intitle:review)+AND+(intext:automatic+OR+intext:automated+OR+intext:automating)+AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)&amp;hl=fr&amp;num=20&amp;as_sdt=0,5" class="gs_or_nvi">En cache</a> <a href="javascript:void(0)" title="Moins" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_DMP21

def test_article_B11():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="M_CxnKh_7awJ" data-did="M_CxnKh_7awJ" data-lid="" data-aid="M_CxnKh_7awJ" data-rp="45">
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><span class="gs_ctu"><span class="gs_ct1">[CITATION]</span><span class="gs_ct2">[C]</span></span> <span id="M_CxnKh_7awJ"><b>Automatic Assessment </b>of <b>Java </b>Programming Patterns for Novices</span></h3>
                <div class="gs_a">A Bagini&nbsp;- University of Western Australia, 2011</div>
                <div class="gs_rs"></div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Cite</span></a> <a href="https://scholar.google.com/scholar?cites=12460756106164170803&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=en">Cited by 3</a> <a href="https://scholar.google.com/scholar?q=related:M_CxnKh_7awJ:scholar.google.com/&amp;scioq=+(intext:automatic+OR+intext:automated+OR+intext:automating)+AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:feedback+OR+intitle:CBA)+AND+(intitle:java)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)&amp;hl=en&amp;as_sdt=0,5">Related articles</a> <a href="javascript:void(0)" title="More" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_B11

def test_article_SE12():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="TWov4KHkqCoJ" data-did="TWov4KHkqCoJ" data-lid="" data-aid="TWov4KHkqCoJ" data-rp="21">
              <div class="gs_ggs gs_fl">
                <div class="gs_ggsd">
                  <div class="gs_or_ggsm" ontouchstart="gs_evt_dsp(event)" tabindex="-1"><a href="https://www.scirp.org/html/6-9601091_17557.htm" data-clk="hl=en&amp;sa=T&amp;oi=gga&amp;ct=gga&amp;cd=21&amp;d=3073958129582434893&amp;ei=72CYZZv7OJesy9YPv4Ok2Ao" data-clk-atid="TWov4KHkqCoJ"><span class="gs_ctg2">[HTML]</span> scirp.org</a></div>
                </div>
              </div>
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><span class="gs_ctu"><span class="gs_ct1">[CITATION]</span><span class="gs_ct2">[C]</span></span> <span id="TWov4KHkqCoJ">An intelligent <b>assessment </b>tool for students' <b>Java </b>submissions in introductory programming courses</span></h3>
                <div class="gs_a">FA Shamsi, <a href="https://scholar.google.com/citations?user=XU9vJRQAAAAJ&amp;hl=en&amp;oi=sra">A Elnagar</a>&nbsp;- Journal of Intelligent learning systems and applications, 2012</div>
                <div class="gs_rs"></div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Cite</span></a> <a href="https://scholar.google.com/scholar?cites=3073958129582434893&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=en">Cited by 44</a> <a href="https://scholar.google.com/scholar?q=related:TWov4KHkqCoJ:scholar.google.com/&amp;scioq=+(intext:automatic+OR+intext:automated+OR+intext:automating)+AND+(intitle:grading+OR+intitle:evaluation+OR+intitle:correction+OR+intitle:assessment+OR+intitle:feedback+OR+intitle:CBA)+AND+(intitle:java)+AND+(intext:library+OR+intext:open-source+OR+intext:%22open+source%22+OR+intext:libre)+AND+(intext:student+OR+intext:learn+OR+intext:exam+OR+intext:essay+OR+intext:assignment)&amp;hl=en&amp;as_sdt=0,5">Related articles</a> <a href="https://scholar.google.com/scholar?cluster=3073958129582434893&amp;hl=en&amp;as_sdt=0,5" class="gs_nph">All 10 versions</a> <a href="javascript:void(0)" title="More" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_SE12

def test_article_SGK():
    article_entry = """<div class="gs_r gs_or gs_scl" data-cid="_Jz8oW0mcjkJ" data-did="_Jz8oW0mcjkJ" data-lid="" data-aid="_Jz8oW0mcjkJ" data-rp="116">
              <div class="gs_ri">
                <h3 class="gs_rt" ontouchstart="gs_evt_dsp(event)"><a id="_Jz8oW0mcjkJ" href="https://onlinelibrary.wiley.com/doi/abs/10.1002/smr.2643" data-clk="hl=en&amp;sa=T&amp;ct=res&amp;cd=116&amp;d=4139413259817884924&amp;ei=0nmYZcvCMLvNy9YP_9mm4Ao" data-clk-atid="_Jz8oW0mcjkJ">The vital role of community in <b>open source </b>software development: A framework for <b>assessment </b>and ranking</a></h3>
                <div class="gs_a">J Singh, A Gupta, P Kanwal&nbsp;- Journal of Software: Evolution and&nbsp;… - Wiley Online <b>Library</b></div>
                <div class="gs_rs">… teams, identifies core FOSS <b>evaluation</b> criteria and introduces a fully <b>automated</b> method for <br>
                  … primary data for <b>evaluation</b>. To demonstrate the effectiveness of <b>Open</b> <b>Source</b> Quality Radar, …</div>
                <div class="gs_fl gs_flb"><a href="javascript:void(0)" class="gs_or_sav gs_or_btn" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.5 11.57l3.824 2.308-1.015-4.35 3.379-2.926-4.45-.378L7.5 2.122 5.761 6.224l-4.449.378 3.379 2.926-1.015 4.35z"></path>
                    </svg><span class="gs_or_btn_lbl">Save</span></a> <a href="javascript:void(0)" class="gs_or_cit gs_or_btn gs_nph" role="button" aria-controls="gs_cit" aria-haspopup="true"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M6.5 3.5H1.5V8.5H3.75L1.75 12.5H4.75L6.5 9V3.5zM13.5 3.5H8.5V8.5H10.75L8.75 12.5H11.75L13.5 9V3.5z"></path>
                    </svg><span>Cite</span></a> <a href="javascript:void(0)" title="More" class="gs_or_mor gs_oph" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M0.75 5.5l2-2L7.25 8l-4.5 4.5-2-2L3.25 8zM7.75 5.5l2-2L14.25 8l-4.5 4.5-2-2L10.25 8z"></path>
                    </svg></a> <a href="javascript:void(0)" title="Fewer" class="gs_or_nvi gs_or_mor" role="button"><svg viewBox="0 0 15 16" class="gs_or_svg">
                      <path d="M7.25 5.5l-2-2L0.75 8l4.5 4.5 2-2L4.75 8zM14.25 5.5l-2-2L7.75 8l4.5 4.5 2-2L11.75 8z"></path>
                    </svg></a></div>
              </div>
            </div>"""
    article = Parser.parse_article_entry(article_entry)
    assert article == article_SGK

def test_html():
    h = Path("tests/resources/res1.html")
    articles = Parser((h, )).parse()
    assert len(articles) == 10
    assert article_AAA in articles
    assert article_IAK in articles
